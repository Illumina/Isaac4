/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2017 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file BinningFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <cerrno>
#include <fstream>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "alignment/BinMetadata.hh"
#include "alignment/matchSelector/BinningFragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{


/**
 * \brief creates bin data files at regular genomic intervals.
 *
 * \param contigsPerBinMax maximum number of different contigs to be associated with the
 *                         same data file path
 */
static boost::filesystem::path makeDataFilePath(
    std::size_t genomicOffset,
    const uint64_t targetBinLength,
    unsigned binIndex,
    const reference::ReferencePosition& binStartPos,
    const unsigned contigsPerBinMax,
    const reference::SortedReferenceMetadata::Contigs& contigs,
    const bfs::path& binDirectory,
    unsigned &lastBinLastContig,
    uint64_t& lastBinRoundedGenomicOffset,
    unsigned & lastBinContigs,
    uint64_t& lastFileNameGenomicOffset)
{
    reference::ReferencePosition fileNameReferencePosition;
    uint64_t roundedGenomicOffset = (genomicOffset / targetBinLength) * targetBinLength;
    uint64_t fileNameGenomicOffset = roundedGenomicOffset;

    if (binIndex && lastBinRoundedGenomicOffset == roundedGenomicOffset)
    {
        if (lastBinLastContig != binStartPos.getContigId())
        {
            ++lastBinContigs;
        }

        if (contigsPerBinMax < lastBinContigs)
        {
            fileNameGenomicOffset = genomicOffset;
            lastFileNameGenomicOffset = genomicOffset;
            lastBinContigs = 0;
        }
        else
        {
            fileNameGenomicOffset = lastFileNameGenomicOffset;
        }
    }
    else
    {
        lastBinRoundedGenomicOffset = roundedGenomicOffset;
        lastFileNameGenomicOffset = fileNameGenomicOffset;
        lastBinContigs = 0;
    }
    lastBinLastContig = binIndex ? binStartPos.getContigId() : 0;
    fileNameReferencePosition = reference::genomicOffsetToPosition(fileNameGenomicOffset, contigs);

    //            ISAAC_THREAD_CERR << "roundedGenomicOffset:" << roundedGenomicOffset << std::endl;
    //            ISAAC_THREAD_CERR << "fileNameReferencePosition:" << fileNameReferencePosition << std::endl;
    boost::filesystem::path binPath =
        // Pad file names well, so that we don't have to worry about them becoming of different length.
        // This is important for memory reservation to be stable
        binDirectory / (boost::format("bin-%08d-%09d.dat")
            % (binIndex ? (fileNameReferencePosition.getContigId() + 1) : 0)
            % (binIndex ? fileNameReferencePosition.getPosition() : 0)).str();
    return binPath;
}

static void buildBinPathList(
    const alignment::matchSelector::BinIndexMap &binIndexMap,
    const bfs::path &binDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const reference::SortedReferenceMetadata::Contigs &contigs,
    const uint64_t targetBinLength,
    alignment::BinMetadataList &ret)
{
    ISAAC_TRACE_STAT("before buildBinPathList");
    ISAAC_ASSERT_MSG(!binIndexMap.empty(), "Empty binIndexMap is illegal");
    ISAAC_ASSERT_MSG(!binIndexMap.back().empty(), "Empty binIndexMap entry is illegal" << binIndexMap);

    reference::SortedReferenceMetadata::Contigs offsetOderedContigs = contigs;
    std::sort(offsetOderedContigs.begin(), offsetOderedContigs.end(),
            [](const reference::SortedReferenceMetadata::Contig &left,
               const reference::SortedReferenceMetadata::Contig &right)
               {return left.genomicPosition_ < right.genomicPosition_;});

    std::size_t genomicOffset = 0;
    static const unsigned CONTIGS_PER_BIN_MAX = 16;
    unsigned lastBinContigs = 0;
    unsigned lastBinLastContig = -1;
    uint64_t lastBinRoundedGenomicOffset = uint64_t(0) - 1;
    uint64_t lastFileNameGenomicOffset = uint64_t(0) - 1;
    BOOST_FOREACH(const std::vector<unsigned> &contigBins, binIndexMap)
    {
        ISAAC_ASSERT_MSG(!contigBins.empty(), "Unexpected empty contigBins");
        // this offset goes in bin length increments. They don't sum up to the real contig length
        std::size_t binGenomicOffset = genomicOffset;
        for (unsigned i = contigBins.front(); contigBins.back() >= i; ++i)
        {
            const reference::ReferencePosition binStartPos = binIndexMap.getBinFirstPos(i);
            // binIndexMap contig 0 is unaligned bin
            ISAAC_ASSERT_MSG(!i || binIndexMap.getBinIndex(binStartPos) == i, "BinIndexMap is broken");

            const boost::filesystem::path binPath =
                makeDataFilePath(
                    binGenomicOffset, targetBinLength, i, binStartPos,
                    CONTIGS_PER_BIN_MAX, offsetOderedContigs, binDirectory,
                    lastBinLastContig, lastBinRoundedGenomicOffset,
                    lastBinContigs, lastFileNameGenomicOffset);

            const uint64_t binLength = i ? binIndexMap.getBinFirstInvalidPos(i) - binStartPos : 0;
            ret.push_back(
                alignment::BinMetadata(barcodeMetadataList.size(), ret.size(), binStartPos, binLength, binPath));

//            ISAAC_THREAD_CERR << "binPathList.back():" << binPathList.back() << std::endl;
            binGenomicOffset += binLength;
        }
//        ISAAC_THREAD_CERR << "contigBins.front():" << contigBins.front() << std::endl;
        if (!ret.back().isUnalignedBin())
        {
            genomicOffset += contigs.at(ret.back().getBinStart().getContigId()).totalBases_;
        }
    }
}

BinningFragmentStorage::BinningFragmentStorage(
    const boost::filesystem::path &tempDirectory,
    const bool keepUnaligned,
    const BinIndexMap &binIndexMap,
    const reference::SortedReferenceMetadata::Contigs& contigs,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const bool preAllocateBins,
    const uint64_t expectedBinSize,
    const uint64_t targetBinLength,
    const unsigned threads,
    alignment::BinMetadataList &binMetadataList):
        FragmentBinner(keepUnaligned, binIndexMap, preAllocateBins ? expectedBinSize : 0, threads),
        binIndexMap_(binIndexMap),
        expectedBinSize_(expectedBinSize),
        binMetadataList_(binMetadataList)
{
    buildBinPathList(
        binIndexMap, tempDirectory, barcodeMetadataList, contigs, targetBinLength, binMetadataList_);

    const std::size_t worstCaseEstimatedUnalignedBins =
        // This assumes that none of the data will align and we will need to
        // make as many unaligned bins as we expect the aligned ones to be there
        // This is both bad and weak. Bad because we allocate pile of memory
        // that will never be used with proper aligning data.
        // Weak because if in fact nothing will align, we might need more unaligned bins than we expect.
        // Keeping this here because we're talking about a few thousands relatively small structures,
        // so pile is not large enough to worry about, and the reallocation will occur while no other threads
        // are using the list, so it should not invalidate any references.
        binMetadataList_.size() * binIndexMap.getBinLength() / targetBinLength +
        // in case the above math returns 0, we'll have room for at least one unaligned bin
        1;

    binMetadataList_.reserve(binMetadataList_.size() + worstCaseEstimatedUnalignedBins);

    unalignedBinMetadataReserve_.resize(worstCaseEstimatedUnalignedBins, binMetadataList_.back());

    FragmentBinner::open(binMetadataList_.begin(), binMetadataList_.end());
}

BinningFragmentStorage::~BinningFragmentStorage()
{
}

void BinningFragmentStorage::store(
    const BamTemplate &bamTemplate,
    const unsigned barcodeIdx,
    const unsigned threadNumber)
{
    common::StaticVector<char, READS_MAX * (sizeof(io::FragmentHeader) + FRAGMENT_BYTES_MAX)> buffer;
    if (2 == bamTemplate.getFragmentCount())
    {
        packPairedFragment(bamTemplate, 0, barcodeIdx, binIndexMap_, std::back_inserter(buffer));
        const io::FragmentAccessor &fragment0 = reinterpret_cast<const io::FragmentAccessor&>(buffer.front());

        packPairedFragment(bamTemplate, 1, barcodeIdx, binIndexMap_, std::back_inserter(buffer));
        const io::FragmentAccessor &fragment1 = *reinterpret_cast<const io::FragmentAccessor*>(&buffer.front() + fragment0.getTotalLength());

        storePaired(fragment0, fragment1, binMetadataList_, threadNumber);
    }
    else
    {
        packSingleFragment(bamTemplate, barcodeIdx, std::back_inserter(buffer));
        const io::FragmentAccessor &fragment = reinterpret_cast<const io::FragmentAccessor&>(buffer.front());
        storeSingle(fragment, binMetadataList_, threadNumber);
    }
}

void BinningFragmentStorage::prepareFlush() noexcept
{
    if (binMetadataList_.front().getDataSize() > expectedBinSize_)
    {
        ISAAC_ASSERT_MSG(!unalignedBinMetadataReserve_.empty(), "Unexpectedly ran out of reserved BinMetadata when extending the unaligned bin");
        // does not cause memory allocation because of reserve in constructor
        binMetadataList_.resize(binMetadataList_.size() + 1);
        using std::swap;
        swap(binMetadataList_.back(), unalignedBinMetadataReserve_.back());
        unalignedBinMetadataReserve_.pop_back();
        binMetadataList_.back() = binMetadataList_.front();
        binMetadataList_.front().startNew();
    }
}


} //namespace matchSelector
} // namespace alignment
} // namespace isaac
