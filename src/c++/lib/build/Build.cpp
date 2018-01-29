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
 ** \file Build.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 **
 ** \author Roman Petrovski
 **/

#include "common/config.h"

#ifdef HAVE_NUMA
#include <numa.h>
#include <numaif.h>
#endif //HAVE_NUMA
 
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/regex.hpp>

#include "bam/Bam.hh"
#include "bam/BamIndexer.hh"
#include "bgzf/BgzfCompressor.hh"
#include "build/Build.hh"
#include "build/IndelLoader.hh"
#include "common/Debug.hh"
#include "common/FileSystem.hh"
#include "common/Threads.hpp"
#include "io/Fragment.hh"
#include "reference/ContigLoader.hh"

#include "BuildStatsXml.hh"
#include "SortedReferenceXmlBamHeaderAdapter.hh"

namespace isaac
{
namespace build
{

const unsigned BuildContigMap::UNMAPPED_CONTIG;
/**
 * \return Returns the total memory in bytes required to load the bin data and indexes
 */
static uint64_t getBinTotalSize(const alignment::BinMetadata & binMetadata)
{
    return
        binMetadata.getDataSize() +
        binMetadata.getFIdxElements() * sizeof(FStrandFragmentIndex) +
        binMetadata.getRIdxElements() * sizeof(RStrandOrShadowFragmentIndex) +
        binMetadata.getSeIdxElements() * sizeof(SeFragmentIndex);
}

uint64_t Build::estimateBinCompressedDataRequirements(
    const alignment::BinMetadata & binMetadata,
    const unsigned outputFileIndex) const
{
    // TODO: put the real number in here.
    static const uint64_t EMPTY_BGZF_BLOCK_SIZE = 1234UL;
    if (!binMetadata.getTotalElements())
    {
        return EMPTY_BGZF_BLOCK_SIZE;
    }

    uint64_t thisOutputFileBarcodeElements = 0;
    // accumulate size required to store all barcodes that map to the same output file
    BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
    {
        const unsigned barcodeOutputFileIndex = barcodeBamMapping_.getSampleIndex(barcode.getIndex());
        if (outputFileIndex == barcodeOutputFileIndex)
        {
            const unsigned barcodeIndex = barcode.getIndex();
            ISAAC_ASSERT_MSG(0 != binMetadata.getTotalElements() || 0 == binMetadata.getBarcodeElements(barcodeIndex), "Can't have empty bin with non-empty bin barcode");

            thisOutputFileBarcodeElements += binMetadata.getBarcodeElements(barcodeIndex);
        }
    }

    // assume all data will take the same fraction or less than the number derived from demultiplexed fragments.
    return EMPTY_BGZF_BLOCK_SIZE +
        ((getBinTotalSize(binMetadata) * thisOutputFileBarcodeElements +
            binMetadata.getTotalElements() - 1) / binMetadata.getTotalElements()) * expectedBgzfCompressionRatio_;
}

inline bool orderBySampleIndex(
    const demultiplexing::BarcodePathMap &barcodeBamMapping,
    const flowcell::BarcodeMetadata &left,
    const flowcell::BarcodeMetadata &right)
{
    return barcodeBamMapping.getSampleIndex(left.getIndex()) < barcodeBamMapping.getSampleIndex(right.getIndex());
}

std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > Build::createOutputFileStreams(
    const flowcell::TileMetadataList &tileMetadataList,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    boost::ptr_vector<bam::BamIndex> &bamIndexes) const
{
    unsigned sinkIndexToCreate = 0;
    std::vector<boost::shared_ptr<boost::iostreams::filtering_ostream> > ret;
    ret.reserve(barcodeBamMapping_.getTotalSamples());

    createDirectories(barcodeBamMapping_, barcodeMetadataList);

    flowcell::BarcodeMetadataList barcodesOrderedBySample(barcodeMetadataList);
    std::sort(barcodesOrderedBySample.begin(), barcodesOrderedBySample.end(),
              boost::bind(&orderBySampleIndex, boost::ref(barcodeBamMapping_), _1, _2));
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodesOrderedBySample)
    {
        if (sinkIndexToCreate == barcodeBamMapping_.getSampleIndex(barcode.getIndex()))
        {
            const boost::filesystem::path &bamPath = barcodeBamMapping_.getFilePath(barcode);
            if (!barcode.isUnmappedReference())
            {
                ISAAC_THREAD_CERR << "Created BAM file: " << bamPath << std::endl;

                const reference::SortedReferenceMetadata &sampleReference =
                    sortedReferenceMetadataList_.at(barcode.getReferenceIndex());

                std::string compressedHeader;
                {
                    std::ostringstream oss(compressedHeader);
                    boost::iostreams::filtering_ostream bgzfStream;
                    bgzfStream.push(bgzf::BgzfCompressor(bamGzipLevel_),65535,0);
                    bgzfStream.push(oss);
                    bam::serializeHeader(bgzfStream,
                                         argv_,
                                         description_,
                                         bamHeaderTags_,
                                         bamPuFormat_,
                                         makeSortedReferenceXmlBamHeaderAdapter(
                                             sampleReference,
                                             boost::bind(&BuildContigMap::isMapped, &contigMap_, barcode.getReferenceIndex(), _1),
                                             tileMetadataList, barcodeMetadataList,
                                             barcode.getSampleName()));
                    bgzfStream.strict_sync();
                    compressedHeader = oss.str();
                }

                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>(new boost::iostreams::filtering_ostream()));
                boost::iostreams::filtering_ostream &bamStream = *ret.back();
                if (bamProduceMd5_)
                {
                    bamStream.push(io::FileSinkWithMd5(bamPath.c_str(), std::ios_base::binary));
                }
                else
                {
                    bamStream.push(boost::iostreams::basic_file_sink<char>(bamPath.string(), std::ios_base::binary));
                }

                if (!bamStream) {
                    BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAM file " + bamPath.string()));
                }

                if (!bamStream.write(compressedHeader.c_str(), compressedHeader.size()))
                {
                    BOOST_THROW_EXCEPTION(
                        common::IoException(errno, (boost::format("Failed to write %d bytes into stream %s") %
                            compressedHeader.size() % bamPath.string()).str()));
                }

                // Create BAM Indexer
                unsigned headerCompressedLength = compressedHeader.size();
                unsigned contigCount = sampleReference.getFilteredContigsCount(
                    boost::bind(&BuildContigMap::isMapped, &contigMap_, barcode.getReferenceIndex(), _1));
                bamIndexes.push_back(new bam::BamIndex(bamPath, contigCount, headerCompressedLength));
            }
            else
            {
                ret.push_back(boost::shared_ptr<boost::iostreams::filtering_ostream>());
                bamIndexes.push_back(new bam::BamIndex());
                ISAAC_THREAD_CERR << "Skipped BAM file due to unmapped barcode reference: " << bamPath << " " << barcode << std::endl;
            }
            ++sinkIndexToCreate;
        }
    }
    ISAAC_ASSERT_MSG(barcodeBamMapping_.getTotalSamples() == sinkIndexToCreate, "must create all output file sinks");

    return ret;
}

alignment::BinMetadataCRefList filterBins(
    alignment::BinMetadataList& bins,
    const std::string &binRegexString)
{
    alignment::BinMetadataCRefList ret;

    if ("all" == binRegexString)
    {
        std::transform(bins.begin(), bins.end(), std::back_inserter(ret),
                       [](alignment::BinMetadata &bm){return boost::ref(bm);});
    }
    else if ("skip-empty" == binRegexString)
    {
        BOOST_FOREACH(alignment::BinMetadata &binMetadata, bins)
        {
            if (!binMetadata.isEmpty())
            {
                ret.push_back(boost::ref(binMetadata));
            }
        }
    }
    else // use regex to filter bins by name
    {
        std::string regexString(binRegexString);
        std::replace(regexString.begin(), regexString.end(), ',', '|');
        boost::regex re(regexString);
        BOOST_FOREACH(alignment::BinMetadata &bin, bins)
        {
            if (!bin.isEmpty() && boost::regex_search(bin.getPath().filename().string(), re))
            {
                ret.push_back(boost::ref(bin));
            }
        }
        if (ret.empty())
        {
            ISAAC_THREAD_CERR << "WARNING: Bam files will be empty. No bins are left after applying the following regex filter: "
                << regexString << std::endl;
        }
    }
    return ret;
}

template <bool unalignedFirst>
bool binLess(const alignment::BinMetadata &left, const alignment::BinMetadata &right)
{
    if (left.isUnalignedBin())
    {
        if (right.isUnalignedBin())
        {
            return left.getDataOffset() < right.getDataOffset();
        }
        return unalignedFirst;
    }
    else
    {
        if (right.isUnalignedBin())
        {
            return !unalignedFirst;
        }

        return left.getBinStart() < right.getBinStart();
    }
}

/**
 * \brief Moves unaligned to the front, to the back or erases them
 */
static alignment::BinMetadataCRefList rearrangeBins(
    alignment::BinMetadataCRefList bins,
    const bool keepUnaligned,
    const bool putUnalignedInTheBack)
{
    bins.erase(std::remove_if(
        bins.begin(), bins.end(),
        [](const alignment::BinMetadata &bm){return bm.isEmpty();}), bins.end());

    if (!keepUnaligned)
    {
        bins.erase(std::remove_if(
            bins.begin(), bins.end(),
            [](const alignment::BinMetadata &bm){return bm.isUnalignedBin();}), bins.end());
    }

    if (putUnalignedInTheBack)
    {
        std::sort(bins.begin(), bins.end(), &binLess<false>);
    }
    else
    {
        std::sort(bins.begin(), bins.end(), &binLess<true>);
    }

    return bins;
}

Build::Build(const std::vector<std::string> &argv,
             const std::string &description,
             const flowcell::FlowcellLayoutList &flowcellLayoutList,
             const flowcell::TileMetadataList &tileMetadataList,
             const flowcell::BarcodeMetadataList &barcodeMetadataList,
             const alignment::BinMetadataList &bins,
             const reference::ReferenceMetadataList &referenceMetadataList,
             const std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
             const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
             const reference::ContigLists &contigLists,
             const boost::filesystem::path outputDirectory,
             const unsigned maxLoaders,
             const unsigned maxComputers,
             const unsigned maxSavers,
             const build::GapRealignerMode realignGaps,
             const unsigned realignMapqMin,
             const boost::filesystem::path &knownIndelsPath,
             const int bamGzipLevel,
             const std::string &bamPuFormat,
             const bool bamProduceMd5,
             const std::vector<std::string> &bamHeaderTags,
             const unsigned expectedCoverage,
             const uint64_t targetBinSize,
             const double expectedBgzfCompressionRatio,
             const bool singleLibrarySamples,
             const bool keepDuplicates,
             const bool markDuplicates,
             const bool anchorMate,
             const bool realignGapsVigorously,
             const bool realignDodgyFragments,
             const unsigned realignedGapsPerFragment,
             const bool clipSemialigned,
             const alignment::AlignmentCfg &alignmentCfg,
             const bool loadAllContigs,
             const std::string &binRegexString,
             const unsigned char forcedDodgyAlignmentScore,
             const bool keepUnaligned,
             const bool putUnalignedInTheBack,
             const IncludeTags includeTags,
             const bool pessimisticMapQ)
    :argv_(argv),
     description_(description),
     flowcellLayoutList_(flowcellLayoutList),
     tileMetadataList_(tileMetadataList),
     barcodeMetadataList_(barcodeMetadataList),
     bins_(bins),
     binRefs_(rearrangeBins(
         filterBins(bins_, binRegexString), keepUnaligned, putUnalignedInTheBack)),
     sortedReferenceMetadataList_(sortedReferenceMetadataList),
     contigMap_(barcodeMetadataList_, binRefs_, sortedReferenceMetadataList_, false),//!loadAllContigs && "skip-empty" == binRegexString),
     outputDirectory_(outputDirectory),
     maxLoaders_(maxLoaders),
     maxComputers_(maxComputers),
     maxRealigners_(1),
     allocatedBins_(0),
     maxSavers_(maxSavers),
     bamGzipLevel_(bamGzipLevel),
     bamPuFormat_(bamPuFormat),
     bamProduceMd5_(bamProduceMd5),
     bamHeaderTags_(bamHeaderTags),
     forcedDodgyAlignmentScore_(forcedDodgyAlignmentScore),
     singleLibrarySamples_(singleLibrarySamples),
     keepDuplicates_(keepDuplicates),
     markDuplicates_(markDuplicates),
     anchorMate_(anchorMate),
     realignGapsVigorously_(realignGapsVigorously),
     realignDodgyFragments_(realignDodgyFragments),
     realignedGapsPerFragment_(realignedGapsPerFragment),
     clipSemialigned_(clipSemialigned),
     alignmentCfg_(alignmentCfg),
     realignGaps_(realignGaps),
     realignMapqMin_(realignMapqMin),
     expectedCoverage_(expectedCoverage),
     targetBinSize_(targetBinSize),
     expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio),
     maxReadLength_(getMaxReadLength(flowcellLayoutList_)),
     includeTags_(includeTags),
     pessimisticMapQ_(pessimisticMapQ),
     forceTermination_(false),
     threads_(maxComputers_ + maxLoaders_ + maxSavers_),
     contigLists_(contigLists),
     barcodeBamMapping_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList_, "sorted.bam")),
     bamIndexes_(),
     bamFileStreams_(createOutputFileStreams(tileMetadataList_, barcodeMetadataList_, bamIndexes_)),
     stats_(binRefs_, barcodeMetadataList_),
     threadBgzfBuffers_(threads_.size(), BgzfBuffers(bamFileStreams_.size())),
     threadBgzfStreams_(threads_.size()),
     threadBamIndexParts_(threads_.size()),
     knownIndels_((build::GapRealignerMode::REALIGN_NONE == realignGaps_ || knownIndelsPath.empty()) ?
         gapRealigner::Gaps() : loadIndels(knownIndelsPath, sortedReferenceMetadataList_)),
     gapRealigner_(threads_.size(),
         realignGapsVigorously, realignDodgyFragments, realignedGapsPerFragment, clipSemialigned,
//         alignmentCfg_.normalizedMismatchScore_,
//         alignmentCfg_.normalizedGapOpenScore_,
//         alignmentCfg_.normalizedGapExtendScore_,
//         alignmentCfg_.normalizedMaxGapExtendScore_,
         barcodeMetadataList, barcodeTemplateLengthStatistics, contigLists_),
     binSorter_(singleLibrarySamples_, keepDuplicates_, markDuplicates_, anchorMate_,
               barcodeBamMapping_, barcodeMetadataList_, contigLists_, alignmentCfg_.splitGapLength_)
{
    computeSlotWaitingBins_.reserve(threads_.size());
    while(threadBgzfStreams_.size() < threads_.size())
    {
        threadBgzfStreams_.push_back(new boost::ptr_vector<boost::iostreams::filtering_ostream>(bamFileStreams_.size()));
    }
    while(threadBamIndexParts_.size() < threads_.size())
    {
        threadBamIndexParts_.push_back(new boost::ptr_vector<bam::BamIndexPart>(bamFileStreams_.size()));
    }

    threads_.execute(boost::bind(&Build::allocateThreadData, this, _1));

    // when number of bins is smaller than number of threads, some complete tasks don't get erased before other tasks are added.
    tasks_.reserve(std::max(binRefs_.size(), threads_.size()));

//    testBinsFitInRam();
}

void Build::allocateThreadData(const std::size_t threadNumber)
{
//#ifdef HAVE_NUMA
//    if (common::isNumaAvailable())
//    {
//        uint64_t nodemask = 1UL << common::ThreadVector::getThreadNumaNode();
//        ISAAC_ASSERT_MSG(-1 != set_mempolicy(MPOL_BIND/*|MPOL_F_STATIC_NODES*/, &nodemask, sizeof(nodemask) * 8),
//                         "set_mempolicy for nodemask: " << nodemask <<
//                         " failed, errno: " << errno << ":" << strerror(errno));
//    }
//#endif //HAVE_NUMA
}

void Build::run(common::ScopedMallocBlock &mallocBlock)
{
    alignment::BinMetadataCRefList::iterator nextUnprocessedBinIt(binRefs_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnallocatedBinIt(binRefs_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnloadedBinIt(binRefs_.begin());
    alignment::BinMetadataCRefList::const_iterator nextUnsavedBinIt(binRefs_.begin());

    threads_.execute(boost::bind(&Build::sortBinParallel, this,
                                boost::ref(nextUnprocessedBinIt),
                                boost::ref(nextUnallocatedBinIt),
                                boost::ref(nextUnloadedBinIt),
                                boost::ref(nextUnsavedBinIt),
                                boost::ref(mallocBlock),
                                _1));

    unsigned fileIndex = 0;
    BOOST_FOREACH(const boost::filesystem::path &bamFilePath, barcodeBamMapping_.getPaths())
    {
        // some of the streams are null_sink (that's when reference is unmapped for the sample).
        // this is the simplest way to ignore them...
        std::ostream *stm = bamFileStreams_.at(fileIndex).get();
        if (stm)
        {
            bam::serializeBgzfFooter(*stm);
            stm->flush();
            ISAAC_THREAD_CERR << "BAM file generated: " << bamFilePath.c_str() << "\n";
            bamIndexes_.at(fileIndex).flush();
            ISAAC_THREAD_CERR << "BAM index generated for " << bamFilePath.c_str() << "\n";
        }
        ++fileIndex;
    }
}

void Build::dumpStats(const boost::filesystem::path &statsXmlPath)
{
    BuildStatsXml statsXml(sortedReferenceMetadataList_, binRefs_, barcodeMetadataList_, stats_);
    std::ofstream os(statsXmlPath.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + statsXmlPath.string()));
    }
    statsXml.serialize(os);
}

uint64_t Build::estimateOptimumFragmentsPerBin(
    const unsigned int estimatedFragmentSize,
    const uint64_t availableMemory,
    const double expectedBgzfCompressionRatio,
    const unsigned computeThreads)
{
//    const size_t maxFragmentIndexBytes = std::max(sizeof(io::RStrandOrShadowFragmentIndex),
//                                                         sizeof(io::FStrandFragmentIndex));

    const std::size_t maxFragmentDedupedIndexBytes = sizeof(PackedFragmentBuffer::Index);
    const std::size_t maxFragmentCompressedBytes = estimatedFragmentSize * expectedBgzfCompressionRatio;

    const std::size_t fragmentMemoryRequirements =
        // assume the initial indexes don't stay in memory for too long //maxFragmentIndexBytes              //index containing duplicates
        + estimatedFragmentSize                 //data
        + maxFragmentDedupedIndexBytes     //deduplicated index
        + maxFragmentCompressedBytes       //bgzf chunk
        ;

    // reasonable amount of bins-in-progress to allow for no-delay input/compute/output overlap
//    const unsigned minOverlap = 3;;
    // try to increase granularity so that the CPU gets efficiently utilized. Gap realigner can only
    // use threads of following bins, let's make sure we have plenty of bins allocated so that threare
    // are always some threads ready to help with realignment
    const unsigned minOverlap = computeThreads * 2;
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin estimatedFragmentSize: " << estimatedFragmentSize << "\n";
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin maxFragmentDedupedIndexBytes: " << maxFragmentDedupedIndexBytes << "\n";
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin maxFragmentCompressedBytes: " << maxFragmentCompressedBytes << "\n";
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin availableMemory: " << availableMemory << "\n";
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin fragmentMemoryRequirements: " << fragmentMemoryRequirements << "\n";
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin minOverlap: " << minOverlap << "\n";
    ISAAC_THREAD_CERR << "estimateOptimumFragmentsPerBin availableMemory / fragmentMemoryRequirements / minOverlap: " << (availableMemory / fragmentMemoryRequirements / minOverlap) << "\n";
    return availableMemory / fragmentMemoryRequirements / minOverlap;
}

/**
 * \brief Attempts to reserve memory buffers required to process a bin.
 *
 * \return Non-zero size in bytes if the reservation failed.
 */
void Build::reserveBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadata &bin,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber,
    boost::shared_ptr<BinData> &binDataPtr)
{
    common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams = threadBgzfStreams_.at(threadNumber);
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts = threadBamIndexParts_.at(threadNumber);
    // bin stats have an entry per filtered bin reference.
    const unsigned binStatsIndex = std::distance<alignment::BinMetadataCRefList::const_iterator>(binRefs_.begin(), thisThreadBinIt);
    common::ScopedMallocBlockUnblock unblockMalloc(mallocBlock);
    ISAAC_TRACE_STAT("Before allocating data for " << bin);
    reserveBuffers(
        bin, binStatsIndex, contigLists_, bgzfStreams, bamIndexParts,
        threadBgzfBuffers_.at(threadNumber), binDataPtr);
    ISAAC_TRACE_STAT("After  allocating data for " << bin);
}

void Build::cleanupBinAllocationFailure(
    const alignment::BinMetadata& bin,
    boost::ptr_vector<boost::iostreams::filtering_ostream>& bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart>& bamIndexParts,
    boost::shared_ptr<BinData>& binDataPtr, BgzfBuffers& bgzfBuffers)
{
    bgzfStreams.clear();
    bamIndexParts.clear();
    // give a chance other threads to allocate what they need... TODO: this is not required anymore as allocation happens orderly
    binDataPtr.reset();
    for(bam::BgzfBuffer &bgzfBuffer : bgzfBuffers)
    {
        bam::BgzfBuffer().swap(bgzfBuffer);
    }
    // reset errno, to prevent misleading error messages when failing code does not set errno
    errno = 0;
}

/**
 * \brief Attempts to reserve memory buffers required to process a bin.
 *
 * \return Non-zero size in bytes if the reservation failed.
 */
void Build::reserveBuffers(
    const alignment::BinMetadata &bin,
    const unsigned binStatsIndex,
    const reference::ContigLists &contigLists,
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts,
    BgzfBuffers &bgzfBuffers,
    boost::shared_ptr<BinData> &binDataPtr)
{
    try
    {
        binDataPtr = boost::shared_ptr<BinData>(
            new BinData(realignedGapsPerFragment_,
                        barcodeBamMapping_, barcodeMetadataList_,
                        realignGaps_, realignMapqMin_, knownIndels_, bin, binStatsIndex, tileMetadataList_, contigMap_, contigLists, maxReadLength_,
                        forcedDodgyAlignmentScore_,  flowcellLayoutList_, includeTags_, pessimisticMapQ_, alignmentCfg_.splitGapLength_,
                        expectedCoverage_));

        unsigned outputFileIndex = 0;
        for(bam::BgzfBuffer &bgzfBuffer : bgzfBuffers)
        {
            bgzfBuffer.reserve(estimateBinCompressedDataRequirements(bin, outputFileIndex++));
        }

        ISAAC_ASSERT_MSG(!bgzfStreams.size(), "Expecting empty pool of streams");
        while(bgzfStreams.size() < bamFileStreams_.size())
        {
            bgzfStreams.push_back(new boost::iostreams::filtering_ostream);
            bgzfStreams.back().push(bgzf::BgzfCompressor(bamGzipLevel_), 65535, 0);
            bgzfStreams.back().push(
                boost::iostreams::back_insert_device<bam::BgzfBuffer >(
                    bgzfBuffers.at(bgzfStreams.size()-1)));
            bgzfStreams.back().exceptions(std::ios_base::badbit);
        }

        ISAAC_ASSERT_MSG(!bamIndexParts.size(), "Expecting empty pool of bam index parts");
        while(bamIndexParts.size() < bamFileStreams_.size())
        {
            bamIndexParts.push_back(new bam::BamIndexPart);
        }
    }
    catch (...)
    {
        cleanupBinAllocationFailure(bin, bgzfStreams, bamIndexParts, binDataPtr, bgzfBuffers);
        throw;
    }
}

template <typename ExceptionType>
const char *getExceptionName(ExceptionType &e)
{
    try
    {
        throw e;
    }
    catch (std::bad_alloc &)
    {
        return "std::bad_alloc";
    }
    catch (boost::iostreams::zlib_error &)
    {
        return "boost::iostreams::zlib_error";
    }
    catch (common::IoException &)
    {
        return "common::IoException";
    }
    catch (...)
    {
        return "Unknown Exception";
    }

}

template <typename ExceptionType, typename ExceptionDataT>
bool Build::handleBinAllocationFailure(
    bool warningTraced,
    const alignment::BinMetadata &bin,
    const ExceptionType &e,
    const ExceptionDataT &errorData)
{
    if (forceTermination_)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }

    if (!allocatedBins_)
    {
        forceTermination_ = true;
        stateChangedCondition_.notify_all();
        // couldn't allocate our bin and not bins currently allocated. No way we will ever be able to allocate it
        BOOST_THROW_EXCEPTION(common::ThreadingException(
            (boost::format("ERROR: Failing due to: %s blocking everything with %s : %s Error data: %s ")
                        % bin % getExceptionName(e) % e.what() % errorData).str()));
    }

    if (!warningTraced)
    {
        ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
            bin << " until " << e.what() << " clears. Error data: " << errorData << std::endl;
        warningTraced = true;
    }

    return warningTraced;
}

boost::shared_ptr<BinData> Build::allocateBin(
    boost::unique_lock<boost::mutex> &lock,
    alignment::BinMetadataCRefList::iterator &thisThreadBinsEndIt,
    alignment::BinMetadataCRefList::iterator &nextUnprocessedBinIt,
    alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
    const alignment::BinMetadataCRefList::const_iterator binsEnd,
    common::ScopedMallocBlock &mallocBlock,
    const std::size_t threadNumber)
{
    bool warningTraced = false;

    alignment::BinMetadataCRefList::iterator thisThreadBinIt = thisThreadBinsEndIt;
    alignment::BinMetadata &bin = *thisThreadBinIt;

    // decide what part of data we can process on this thread
    while (++thisThreadBinsEndIt != binsEnd &&
        thisThreadBinIt->get().sameContig(*thisThreadBinsEndIt) &&
        thisThreadBinIt->get().samePath(*thisThreadBinsEndIt) &&
        targetBinSize_ > bin.getDataSize() + thisThreadBinsEndIt->get().getDataSize())
    {
        bin.merge(*thisThreadBinsEndIt);
    }
    ISAAC_THREAD_CERR << "MERGED targetBinSize_:" << targetBinSize_ << " " << bin <<
        " " << bin.getTotalSplitCount() << std::endl;

    // now next thread can start taking another chunk of bins
    nextUnprocessedBinIt = thisThreadBinsEndIt;

    boost::shared_ptr<BinData> ret(0);
    while(nextUnallocatedBinIt != thisThreadBinIt)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        if (!yieldIfPossible(lock, threadNumber, 0))
        {
            stateChangedCondition_.wait(lock);
        }
    }

    while(true)
    {
        try
        {
            reserveBuffers(lock, thisThreadBinIt, bin, mallocBlock, threadNumber, ret);
            break;
        }
        catch (std::bad_alloc &a)
        {
            uint64_t totalBuffersNeeded = 0UL;
            for(unsigned outputFileIndex = 0; outputFileIndex < threadBgzfStreams_.at(threadNumber).size(); ++outputFileIndex)
            {
                totalBuffersNeeded += estimateBinCompressedDataRequirements(bin, outputFileIndex++);
            }
            warningTraced = handleBinAllocationFailure(
                warningTraced, bin, a, BinData::getMemoryRequirements(bin) + totalBuffersNeeded);
        }
        catch (boost::iostreams::zlib_error &z)
        {
            warningTraced = handleBinAllocationFailure(warningTraced, bin, z, z.error());
        }
        catch (common::IoException &io)
        {
            if (EMFILE == io.getErrorNumber())
            {
                warningTraced = handleBinAllocationFailure(warningTraced, bin, io, "Increase number of open files for better efficiency");
            }
            else
            {
                throw;
            }
        }

        stateChangedCondition_.wait(lock);
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
    }

    ++allocatedBins_;

    nextUnallocatedBinIt = thisThreadBinsEndIt;
    return ret;
}

void Build::waitForLoadSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinsEndIt,
    alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt)
{
    bool warningTraced = false;

    while(nextUnloadedBinIt != thisThreadBinIt || !maxLoaders_)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        if (nextUnloadedBinIt == thisThreadBinIt && !warningTraced)
        {
            ISAAC_THREAD_CERR << "WARNING: Holding up processing of bin: " <<
                thisThreadBinIt->get().getPath().c_str() << " until a load slot is available" << std::endl;
            warningTraced = true;
        }

        stateChangedCondition_.wait(lock);
    }

    nextUnloadedBinIt = thisThreadBinsEndIt;
    --maxLoaders_;
}

void Build::returnLoadSlot(const bool exceptionUnwinding)
{
    ++maxLoaders_;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

/**
 * @return true if this thread was the first to set task to 'complete" state
 */
bool Build::executePreemptTask(
    boost::unique_lock<boost::mutex>& lock,
    Task &task,
    const unsigned threadNumber) const
{
    //    ISAAC_THREAD_CERR << "preempt task " << task << " on thread " << threadNumber << "in:" << task->threadsIn_ << std::endl;
    ++task.threadsIn_;
    //    ISAAC_THREAD_CERR << "preempt " << &lock << " " << threadNumber << std::endl;
    try
    {
        task.execute(lock, threadNumber);
    }
    catch (...)
    {
        --task.threadsIn_;
        throw;
    }

    bool ret = false;

    if (!task.complete_)
    {
        //Threads don't come out of execute until there is nothing left to do.
        // stop new threads entering the task;
        task.complete_ = true;
        ret = true;
    }
    --task.threadsIn_;
    //    ISAAC_THREAD_CERR << "preempt task " << task << " on thread " << threadNumber << "in:" << task->threadsIn_ << " done" << std::endl;
    return ret || !task.threadsIn_;
}

bool Build::processMostUrgent(boost::unique_lock<boost::mutex> &lock, const unsigned threadNumber, Task *ownTask)
{
    // find the lowest priority incomplete and not busy task that is not higher than ownTask unless ownTask is 0
    Tasks::iterator highesttPriorityTask = std::min_element(
        tasks_.begin(), tasks_.end(), [ownTask](const Task *left, const Task *right)
        {
            if (ownTask && left->priority_ > ownTask->priority_)
            {
                return false;
            }

            if (left->complete_ || left->busy())
            {
                return false;
            }

            if (ownTask && right->priority_ > ownTask->priority_)
            {
                return true;
            }

            if (right->complete_ || right->busy())
            {
                return true;
            }
            return left->priority_ < right->priority_;
        });

    if(tasks_.end() == highesttPriorityTask)
    {
        return false;
    }

    Task *task = *highesttPriorityTask;
    if(task->complete_ || task->busy() || (ownTask && ownTask->busy() && task->priority_ > ownTask->priority_))
    {
        return false;
    }

    ISAAC_ASSERT_MSG(!ownTask || task->priority_ <= ownTask->priority_, "invalid task found");

    ISAAC_ASSERT_MSG(maxComputers_, "Unexpected maxComputers_ 0");
    --maxComputers_;

    bool ret = false;
    ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnComputeSlot, this, _1))
    {
        //    ISAAC_THREAD_CERR << "preempt task " << task << " on thread " << threadNumber << "in:" << task->threadsIn_ << std::endl;
        ret = executePreemptTask(lock, *task, threadNumber);
    }

    if (ret)
    {
        stateChangedCondition_.notify_all();
    }

    return ret;
}

bool Build::yieldIfPossible(
    boost::unique_lock<boost::mutex>& lock,
    const std::size_t threadNumber,
    Task *task)
{
    bool stateMightHaveChanged = false;
    if (maxComputers_)
    {
        stateMightHaveChanged = processMostUrgent(lock, threadNumber, task);
    }
    return stateMightHaveChanged;
}

template <typename OperationT>
void Build::preemptComputeSlot(
    boost::unique_lock<boost::mutex> &lock,
    const std::size_t maxThreads,
    const std::size_t priority,
    OperationT operation,
    const unsigned threadNumber)
{
    struct OperationTask : public Task
    {
        OperationT operation_;
        OperationTask(const std::size_t maxThreads, const std::size_t priority, OperationT operation):
            Task(maxThreads, priority), operation_(operation) {}
        virtual void execute(boost::unique_lock<boost::mutex> &l, const unsigned tn)
        {
//            ISAAC_THREAD_CERR << "Task::execute " << &l << " " << tn << std::endl;
            operation_(l, tn);
        }
    };
    OperationTask ourTask(maxThreads, priority, operation);
    ISAAC_ASSERT_MSG(tasks_.size() < tasks_.capacity(), "Unexpected high number of concurrent tasks. capacity: " << tasks_.capacity());
    tasks_.push_back(&ourTask);

    stateChangedCondition_.notify_all();

//    ISAAC_THREAD_CERR << "preemptComputeSlot " << &lock << " " << threadNumber << std::endl;
    // keep working until our task is complete.
    while (!ourTask.complete_)
    {
        if (forceTermination_)
        {
            // don't admit new threads
            ourTask.complete_ = true;
        }
        else if (!yieldIfPossible(lock, threadNumber, &ourTask))
        {
            stateChangedCondition_.wait(lock);
        }
    }

    // don't leave while there are some threads still in.
    while (ourTask.threadsIn_)
    {
        stateChangedCondition_.wait(lock);
    }

    ISAAC_ASSERT_MSG(tasks_.end() != std::find(tasks_.begin(), tasks_.end(), &ourTask), "Our task is gone from the list.");
    tasks_.erase(std::find(tasks_.begin(), tasks_.end(), &ourTask));

    if (forceTermination_)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }
}

void Build::returnComputeSlot(const bool exceptionUnwinding)
{
    ++maxComputers_;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::waitForSaveSlot(
    boost::unique_lock<boost::mutex> &lock,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinIt,
    alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt)
{
    while(nextUnsavedBinIt != thisThreadBinIt)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        stateChangedCondition_.wait(lock);
    }
}

void Build::returnSaveSlot(
    alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt,
    const alignment::BinMetadataCRefList::const_iterator thisThreadBinEndIt,
    const bool exceptionUnwinding)
{
    nextUnsavedBinIt = thisThreadBinEndIt;
    if (exceptionUnwinding)
    {
        forceTermination_ = true;
    }
    stateChangedCondition_.notify_all();
}

void Build::sortBinParallel(alignment::BinMetadataCRefList::iterator &nextUnprocessedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnallocatedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnloadedBinIt,
                            alignment::BinMetadataCRefList::const_iterator &nextUnsavedBinIt,
                            common::ScopedMallocBlock &mallocBlock,
                            const std::size_t threadNumber)
{
    boost::unique_lock<boost::mutex> lock(stateMutex_);
    while(binRefs_.end() != nextUnprocessedBinIt)
    {
        static unsigned loadingThreads = 0;
        static unsigned dedupingThreads = 0;
        static unsigned realigningThreads = 0;
        static unsigned serializingThreads = 0;
        static unsigned savingThreads = 0;

        alignment::BinMetadataCRefList::iterator thisThreadBinIt = nextUnprocessedBinIt;
        alignment::BinMetadataCRefList::iterator thisThreadBinsEndIt = thisThreadBinIt;

        // wait and allocate memory required for loading and compressing this bin
        boost::shared_ptr<BinData> binDataPtr =
            allocateBin(lock, thisThreadBinsEndIt, nextUnprocessedBinIt, nextUnallocatedBinIt, binRefs_.end(), mallocBlock, threadNumber);
        waitForLoadSlot(lock, thisThreadBinIt, thisThreadBinsEndIt, nextUnloadedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnLoadSlot, this, _1))
        {
            ++loadingThreads;
    //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                BinLoader binLoader;
                binLoader.loadData(*binDataPtr);
            }
            --loadingThreads;
    //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
        }

        {
            preemptComputeSlot(
                lock, 1, std::distance(binRefs_.begin(), thisThreadBinIt),
                [this, &binDataPtr](boost::unique_lock<boost::mutex> &l, const unsigned tn)
                {
                    ++dedupingThreads;
            //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
                    {
                        common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(l);
                        binSorter_.resolveDuplicates(*binDataPtr, stats_);
                    }
                    --dedupingThreads;
            //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
                },
                threadNumber);

            if (!binDataPtr->isUnalignedBin() && REALIGN_NONE != realignGaps_)
            {
                BinData::iterator nextUnprocessed = binDataPtr->indexBegin();
                int threadsIn = 0;
                preemptComputeSlot(
                    lock, -1, std::distance(binRefs_.begin(), thisThreadBinIt),
                    [this, &threadsIn, &binDataPtr, &nextUnprocessed](boost::unique_lock<boost::mutex> &l, const unsigned tn)
                    {
                        ++threadsIn;
                        ++realigningThreads;
                //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
                        if (nextUnprocessed  == binDataPtr->indexBegin())
                        {
                            ISAAC_THREAD_CERR << "Realigning against " << getTotalGapsCount(binDataPtr->realignerGaps_) <<
                                " unique gaps. " << binDataPtr->bin_ << std::endl;
                        }
                        gapRealigner_.threadRealignGaps(l, *binDataPtr, nextUnprocessed, tn);
                        if (!--threadsIn)
                        {
                            ISAAC_THREAD_CERR << "Realigning gaps done. " << binDataPtr->bin_ << std::endl;
                        }
                        --realigningThreads;
                //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
                    },
                    threadNumber);

            }

            preemptComputeSlot(
                lock, 1, std::distance(binRefs_.begin(), thisThreadBinIt),
                [this, &binDataPtr, &threadNumber](boost::unique_lock<boost::mutex> &l, const unsigned tn)
                {
                    ++serializingThreads;
            //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
                    {
                        common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(l);
                        // Don't use tn!!! the streams have been allocated for the threadNumber.
                        binSorter_.serialize(
                            *binDataPtr, threadBgzfStreams_.at(threadNumber), threadBamIndexParts_.at(threadNumber));
                        threadBgzfStreams_.at(threadNumber).clear();
                    }
                    --serializingThreads;
            //        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
                },
                threadNumber);
        }
        // give back some memory to allow other threads to load
        // data while we're waiting for our turn to save
        binDataPtr.reset();
        stateChangedCondition_.notify_all();

        ++savingThreads;
//        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
        // wait for our turn to store bam data
        waitForSaveSlot(lock, thisThreadBinIt, nextUnsavedBinIt);
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&Build::returnSaveSlot, this, boost::ref(nextUnsavedBinIt), thisThreadBinsEndIt, _1))
        {
            saveAndReleaseBuffers(lock, thisThreadBinIt->get().getPath(), threadNumber);
        }
        --savingThreads;
//        ISAAC_THREAD_CERR << "Threads:" << allocatedBins_ << "," << dedupingThreads << "," << realigningThreads << "," << serializingThreads << "," << savingThreads << "," << loadingThreads << std::endl;
    }

    // Don't release thread until all saving is done. Use threads that don't get anything to process for preemptive tasks such as realignment.
    while(!forceTermination_ && binRefs_.end() != nextUnsavedBinIt)
    {
        if (!yieldIfPossible(lock, threadNumber, 0))
        {
            stateChangedCondition_.wait(lock);
        }
    }
}

/**
 * \brief Save bgzf compressed buffers into corresponding sample files and and release associated memory
 */
void Build::saveAndReleaseBuffers(
    boost::unique_lock<boost::mutex> &lock,
    const boost::filesystem::path &filePath,
    const std::size_t threadNumber)
{
    unsigned index = 0;
    BOOST_FOREACH(bam::BgzfBuffer &bgzfBuffer, threadBgzfBuffers_.at(threadNumber))
    {
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
//            ISAAC_THREAD_CERR << "Freezing thread " << threadNumber << " for bin " << filePath.c_str() << std::endl;
//            while (true)
//            {
//                sleep(1);
//            }
            std::ostream *stm = bamFileStreams_.at(index).get();
            if (!stm)
            {
                ISAAC_ASSERT_MSG(bgzfBuffer.empty(), "Unexpected data for bam file belonging to a sample with unmapped reference");
            }
            else
            {
                saveBuffer(bgzfBuffer, *stm, threadBamIndexParts_.at(threadNumber).at(index), bamIndexes_.at(index), filePath);
            }
        }
        // release rest of the memory that was reserved for this bin
        bam::BgzfBuffer().swap(bgzfBuffer);
        ++index;
    }
    --allocatedBins_;
    threadBamIndexParts_.at(threadNumber).clear();
}

void Build::saveBuffer(
    const bam::BgzfBuffer &bgzfBuffer,
    std::ostream &bamStream,
    const bam::BamIndexPart &bamIndexPart,
    bam::BamIndex &bamIndex,
    const boost::filesystem::path &filePath)
{
    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath.c_str() << std::endl;
    const clock_t start = clock();
    if(!bgzfBuffer.empty() && !bamStream.write(&bgzfBuffer.front(), bgzfBuffer.size())/* ||
        !bamStream.strict_sync()*/){
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed to write bgzf block of %d bytes into bam stream") % bgzfBuffer.size()).str()));
    }
    bamIndex.processIndexPart( bamIndexPart, bgzfBuffer );

    ISAAC_THREAD_CERR << "Saving " << bgzfBuffer.size() << " bytes of sorted data for bin " << filePath.c_str() << " done in " << (clock() - start) / 1000 << "ms\n";
}

} // namespace build
} // namespace isaac
