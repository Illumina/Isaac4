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
 ** \file FragmentBinner.hh
 **
 ** \brief Stores fragments in bin files.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_BINNER_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_BINNER_HH

#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/thread/mutex.hpp>

#include "alignment/BinMetadata.hh"
#include "BinIndexMap.hh"
#include "common/Memory.hh"
#include "io/FileBufCache.hh"
#include "io/Fragment.hh"


namespace isaac
{
namespace alignment
{
namespace matchSelector
{

namespace bfs = boost::filesystem;

class FragmentPacker: boost::noncopyable
{
public:
    template <typename InsertIT>
    static InsertIT packPairedFragment(
        const alignment::BamTemplate &bamTemplate,
        const unsigned fragmentIndex,
        const unsigned barcodeIdx,
        const BinIndexMap &binIndexMap,
        InsertIT insertIt)
    {
        ISAAC_ASSERT_MSG(READS_MAX == bamTemplate.getFragmentCount(), "Expected paired data");

        const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(fragmentIndex);
        const alignment::FragmentMetadata &mate = bamTemplate.getMateFragmentMetadata(fragment);

//        const unsigned mateStorageBin = mate.isNoMatch() && fragment.isNoMatch() ?
//            0 :
//            binIndexMap.getBinIndex(mate.getFStrandReferencePosition());
        // reads are currently stored in every bin that they cover. This means dupe detection will see
        // the reverse alignment duplicate candidates even if they begin in different bins.
        const unsigned mateStorageBin = 0;

        const io::FragmentHeader header(bamTemplate, fragment, mate, barcodeIdx, mateStorageBin);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(header.clusterId_, "FragmentPacker::packPairedFragment " << header << " " << fragment);

        insertIt = std::copy(header.bytesBegin(), header.bytesEnd(), insertIt);
        insertIt = storeBclAndCigar(fragment, insertIt);
        insertIt = std::copy(bamTemplate.nameBegin(), bamTemplate.nameEnd(), insertIt);
        *insertIt++ = 0;

//        ISAAC_ASSERT_MSG("FC:10201664" != std::string(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()), fragment);

        return insertIt;
    }

    template <typename InsertIT>
    static InsertIT packSingleFragment(
        const alignment::BamTemplate &bamTemplate,
        const unsigned barcodeIdx,
        InsertIT insertIt)
    {
        ISAAC_ASSERT_MSG(1 == bamTemplate.getFragmentCount(), "Expected single-ended data");

        const alignment::FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(0);
        const io::FragmentHeader header(bamTemplate, fragment, barcodeIdx);
        insertIt = std::copy(header.bytesBegin(), header.bytesEnd(), insertIt);
        insertIt = storeBclAndCigar(fragment, insertIt);
        insertIt = std::copy(bamTemplate.nameBegin(), bamTemplate.nameEnd(), insertIt);
        *insertIt++ = 0;

        return insertIt;
    }

private:
    static const unsigned READS_MAX = 2;

    template <typename InsertIT>
    static InsertIT storeBclAndCigar(
        const alignment::FragmentMetadata & fragment,
        InsertIT insertIt)
    {
        // copy the bcl data (reverse-complement the sequence if the fragment is reverse-aligned)
        BclClusters::const_iterator bclData = fragment.getBclData();
        if (fragment.isReverse())
        {
            insertIt = std::transform(std::reverse_iterator<BclClusters::const_iterator>(bclData + fragment.getReadLength()),
                                          std::reverse_iterator<BclClusters::const_iterator>(bclData),
                                          insertIt, oligo::getReverseBcl);
        }
        else
        {
            insertIt = std::copy(bclData, bclData + fragment.getReadLength(), insertIt);

        }

        if (fragment.isAligned())
        {
            const alignment::Cigar::const_iterator cigarBegin = fragment.cigarBuffer->begin() + fragment.cigarOffset;
            const alignment::Cigar::const_iterator cigarEnd = cigarBegin + fragment.cigarLength;
            BOOST_FOREACH(const unsigned cig, std::make_pair(cigarBegin, cigarEnd))
            {
                *insertIt++ =cig;
                *insertIt++ =(cig >> 8);
                *insertIt++ =(cig >> 16);
                *insertIt++ =(cig >> 24);
            }
        }
        return insertIt;
    }
};

class FragmentBinner: boost::noncopyable
{
public:
    FragmentBinner(
        const bool keepUnaligned,
        const BinIndexMap &binIndexMap,
        const uint64_t expectedBinSize,
        const unsigned threads);

    // opens a range of bins. Multiple opens are called over the lifetime of FragmentBinner
    void open(
        const alignment::BinMetadataList::iterator binsBegin,
        const alignment::BinMetadataList::iterator binsEnd);

    /**
     * \brief reclaims any unused storage.
     */
    void close() noexcept;

    void storeSingle(
        const io::FragmentAccessor &fragment,
        alignment::BinMetadataList &binMetadataList,
        const unsigned threadNumber);

    void storePaired(
        const io::FragmentAccessor &fragment0,
        const io::FragmentAccessor &fragment1,
        alignment::BinMetadataList &binMetadataList,
        const unsigned threadNumber);

    /**
     * \brief empty all thread buffers and update the corresponding BinMetadata
     */
    void flush(BinMetadataList &binMetadataList);

private:
    /// Maximum number of bins a fragment is expected to cover. In theory this can be up to total number of bins.
    static const unsigned FRAGMENT_BINS_MAX = 10*1024;
    static const unsigned CLUSTER_BINS_MAX = FRAGMENT_BINS_MAX * 2;
    // bytes to store before flushing
    static const std::size_t BUFFER_BYTES_MAX =  4096;
    static const unsigned UNMAPPED_BIN = -1U;

    static const unsigned READS_MAX = 2;
    const bool keepUnaligned_;
    const uint64_t expectedBinSize_;

    const BinIndexMap &binIndexMap_;

    uint64_t binZeroRecordsBinned_;
    // The purpose is to reduce the synchronisation collisions between threads requesting write. So, small number is bad.
    // Large number of mutexes is something to consider. In particular one mutex per file might be the best choice
    // Right now with mutex size of 40 bytes this keeps all of them in one page.
    boost::array<boost::mutex, 4096 / sizeof(boost::mutex)> binMutex_;
    std::vector<io::FileBufWithReopen> files_;
    std::vector<int> binFiles_;

    typedef common::StaticVector<unsigned, CLUSTER_BINS_MAX> FragmentBins;

    // buffer writes on each thread so that the locking and random file writes are a lesser issue
    // guaranteed to fit the bin list, so overflows only on adding new fragments
    typedef common::StaticVector<char, BUFFER_BYTES_MAX + CLUSTER_BINS_MAX * sizeof(unsigned)> FileBuffer;
    typedef std::vector<FileBuffer> FileBuffers;
    std::vector<FileBuffers> threadFileBuffers_;

    static void bufferBinIndexes(
        const FragmentBins &bins,
        FileBuffer &buffer);

    static bool bufferFragment(
        const io::FragmentAccessor &fragment,
        FileBuffer &buffer);

    static bool bufferPair(
        const io::FragmentAccessor &fragment0,
        const io::FragmentAccessor &fragment1,
        FileBuffer &buffer);

    struct BinIndexList
    {
        unsigned indexCount_;
        unsigned indexes_[];
    };

    void flushSingle(
        const io::FragmentAccessor &fragment,
        const BinIndexList &binIndexList,
        alignment::BinMetadataList &binMetadataList,
        const unsigned fileIndex);

    void flushBuffer(
        FileBuffer &buffer,
        alignment::BinMetadataList &binMetadataList,
        const unsigned fileIndex);

    void getFragmentStorageBins(const io::FragmentAccessor &fragment, FragmentBins &bins);

    void reopenBin(const BinMetadata &binMetadata, std::size_t file);
    void openBinFile(const BinMetadata &binMetadata, std::size_t file);
    void registerFragment(const io::FragmentAccessor& fragment,
                          const bool splitRead, const bool realignableSplit, BinMetadata& binMetadata);
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_BINNER_HH
