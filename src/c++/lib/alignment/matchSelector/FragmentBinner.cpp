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
 ** \file FragmentBinner.cpp
 **
 ** \author Roman Petrovski
 **/

#include <cerrno>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>

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

const unsigned FragmentBinner::FRAGMENT_BINS_MAX;
const unsigned FragmentBinner::UNMAPPED_BIN;

FragmentBinner::FragmentBinner(
    const bool keepUnaligned,
    const BinIndexMap &binIndexMap,
    const uint64_t expectedBinSize,
    const unsigned threads):
        keepUnaligned_(keepUnaligned),
        expectedBinSize_(expectedBinSize),
        binIndexMap_(binIndexMap),
        binZeroRecordsBinned_(0),
        threadFileBuffers_(threads)
{
}

void FragmentBinner::registerFragment(const io::FragmentAccessor& fragment,
                                      const bool splitRead,
                                      const bool realignableSplit,
                                      BinMetadata& binMetadata)
{
    if (!fragment.isAligned() && !fragment.isMateAligned())
    {
        // Bin 0 gets split during bam generation. It is important bin 0 chunks reflect distribution in the order
        // in which the records are stored in bin 0. Notice that this is not the case
        // with aligned bins where distribution reflects alignment position
        binMetadata.incrementDataSize(binZeroRecordsBinned_, fragment.getTotalLength());
        binMetadata.incrementNmElements(binZeroRecordsBinned_, 1, fragment.barcode_);
        ++binZeroRecordsBinned_;
    }
    else
    {
        binMetadata.incrementDataSize(fragment.fStrandPosition_, fragment.getTotalLength());
        if (!fragment.flags_.paired_)
        {
            binMetadata.incrementSeIdxElements(fragment.fStrandPosition_, 1, fragment.barcode_);
        }
        else if (fragment.flags_.reverse_ || fragment.flags_.unmapped_)
        {
            binMetadata.incrementRIdxElements(fragment.fStrandPosition_, 1, fragment.barcode_);
        }
        else
        {
            binMetadata.incrementFIdxElements(fragment.fStrandPosition_, 1, fragment.barcode_);
        }

        if (splitRead)
        {
            binMetadata.incrementSplitCount(fragment.fStrandPosition_, realignableSplit, fragment.barcode_);
            binMetadata.incrementGapCount(fragment.fStrandPosition_, realignableSplit, fragment.barcode_);
        }
        else
        {
            binMetadata.incrementGapCount(fragment.fStrandPosition_, fragment.gapCount_, fragment.barcode_);
        }
        binMetadata.incrementCigarLength(
            // approximate aligned bases as readLengt_
            fragment.fStrandPosition_, fragment.cigarLength_, fragment.readLength_, fragment.barcode_);
    }
}

void FragmentBinner::bufferBinIndexes(
    const FragmentBins &bins,
    FileBuffer &buffer)
{
    const typeof(BinIndexList::indexCount_) count = bins.size();
    std::copy(reinterpret_cast<const char *>(&count), reinterpret_cast<const char *>(&count) + sizeof(count), std::back_inserter(buffer));
    std::copy(reinterpret_cast<const char *>(&bins.front()), reinterpret_cast<const char *>(&bins.front()) + bins.size() * sizeof(FragmentBins::value_type), std::back_inserter(buffer));
}

bool FragmentBinner::bufferFragment(
    const io::FragmentAccessor &fragment,
    FileBuffer &buffer)
{
    if (buffer.size() + fragment.getTotalLength() > BUFFER_BYTES_MAX)
    {
        return false;
    }

    std::copy(fragment.begin(), fragment.end(), std::back_inserter(buffer));
    return true;
}

bool FragmentBinner::bufferPair(
    const io::FragmentAccessor &fragment0,
    const io::FragmentAccessor &fragment1,
    FileBuffer &buffer)
{
    const std::size_t before = buffer.size();
    // make sure orphan is stored first in shadow/orphan pair
    if (fragment0.isAligned())
    {
        if (!bufferFragment(fragment0, buffer))
        {
            return false;
        }
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment0.clusterId_, "BinningFragmentStorage::storePaired: " << fragment0);
    }
    if (!bufferFragment(fragment1, buffer))
    {
        buffer.resize(before);
        return false;
    }
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment1.clusterId_, "BinningFragmentStorage::storePaired: " << fragment1);
    if (!fragment0.isAligned())
    {
        if (!bufferFragment(fragment0, buffer))
        {
            buffer.resize(before);
            return false;
        }
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment0.clusterId_, "BinningFragmentStorage::storePaired: " << fragment0);
    }
    return true;
}

void FragmentBinner::flushSingle(
    const io::FragmentAccessor &fragment,
    const BinIndexList &binIndexList,
    alignment::BinMetadataList &binMetadataList,
    const unsigned fileIndex)
{
    unsigned lastBinIndex = -1U;
    for (unsigned i = 0; binIndexList.indexCount_ != i; ++i)
    {
        lastBinIndex = binIndexList.indexes_[i];
        registerFragment(
                fragment,
                // looks like some historical check for unaligned bin. Currently 
                // results in massive undercounting of split alignments. Commented out: //0 != i &&
                fragment.isAligned() && fragment.flags_.splitAlignment_,
                fragment.isAligned() && fragment.flags_.realignableSplit_,
                binMetadataList[lastBinIndex]);
    }

    io::FileBufWithReopen &file = files_.at(fileIndex);

#ifdef ISAAC_TEMP_STORE_DISABLED
    return;
#endif //ISAAC_TEMP_STORE_DISABLED

    if (fragment.getTotalLength() != file.sputn(
            reinterpret_cast<const char*>(&fragment), fragment.getTotalLength()))
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to write into " + binMetadataList[lastBinIndex].getPathString()));
    }
}

void FragmentBinner::flushBuffer(
    FileBuffer &buffer,
    alignment::BinMetadataList &binMetadataList,
    const unsigned fileIndex)
{
//    ISAAC_THREAD_CERR << "flushBuffer fileIndex: " << fileIndex << " for " << buffer.size() << std::endl;
    boost::unique_lock<boost::mutex> lock(binMutex_[fileIndex % binMutex_.size()]);
    for (const char *p = &buffer.front(); &buffer.front() + buffer.size() != p;)
    {
        const io::FragmentAccessor &fragment0 = reinterpret_cast<const io::FragmentAccessor &>(*p);
        ISAAC_ASSERT_MSG(fragment0.flags_.initialized_, "Attempt to store an uninitialised " << fragment0);
        if (fragment0.flags_.paired_)
        {
            const io::FragmentAccessor &fragment1 = *reinterpret_cast<const io::FragmentAccessor *>(fragment0.end());
            ISAAC_ASSERT_MSG(fragment1.flags_.initialized_, "Attempt to store an uninitialised " << fragment1);

            const BinIndexList &binIndexList = *reinterpret_cast<const BinIndexList *>(fragment1.end());
            flushSingle(fragment0, binIndexList, binMetadataList, fileIndex);
            flushSingle(fragment1, binIndexList, binMetadataList, fileIndex);
            p = reinterpret_cast<const char*>(&binIndexList.indexes_[binIndexList.indexCount_]);
        }
        else
        {
            const BinIndexList &binIndexList = *reinterpret_cast<const BinIndexList *>(fragment0.end());
            flushSingle(fragment0, binIndexList, binMetadataList, fileIndex);
            p = reinterpret_cast<const char*>(&binIndexList.indexes_[binIndexList.indexCount_]);
        }
    }
    buffer.clear();
//    ISAAC_THREAD_CERR << "flushBuffer fileIndex: " << fileIndex << " for " << buffer.size() << " done" << std::endl;
}

void FragmentBinner::storePaired(
    const io::FragmentAccessor &fragment0,
    const io::FragmentAccessor &fragment1,
    alignment::BinMetadataList &binMetadataList,
    const unsigned threadNumber)
{
    //return;
    FragmentBins bins;
    getFragmentStorageBins(fragment0, bins);
    getFragmentStorageBins(fragment1, bins);

    FileBuffers &buffers = threadFileBuffers_.at(threadNumber);
    unsigned lastFileIndex = UNMAPPED_BIN;
    FragmentBins lastFileBins;
    for (const unsigned binIndex : bins)
    {
        if (!binIndex && (fragment0.isAligned() || fragment1.isAligned()))
        {
            // only when both reads are unaligned, the pair goes into bin 0
            continue;
        }
        const unsigned fileIndex = binFiles_.at(binIndex);
        if (UNMAPPED_BIN != fileIndex)
        {
            if (lastFileIndex != fileIndex)
            {
                if (UNMAPPED_BIN != lastFileIndex)
                {
                    bufferBinIndexes(lastFileBins, buffers.at(lastFileIndex));
                    lastFileBins.clear();
                }

                FileBuffer &buffer = buffers[fileIndex];
                if (!bufferPair(fragment0, fragment1, buffer))
                {
                    flushBuffer(buffer, binMetadataList, fileIndex);
                    ISAAC_VERIFY_MSG(bufferPair(fragment0, fragment1, buffer), "Could not buffer into empty buffer" << fragment0 << "-" << fragment1);
                }
                lastFileIndex = fileIndex;
            }

            lastFileBins.push_back(binIndex);
        }
    }
    bufferBinIndexes(lastFileBins, buffers.at(lastFileIndex));
}

void FragmentBinner::storeSingle(
    const io::FragmentAccessor &fragment,
    alignment::BinMetadataList &binMetadataList,
    const unsigned threadNumber)
{
    FragmentBins bins;
    getFragmentStorageBins(fragment, bins);

    FileBuffers &buffers = threadFileBuffers_.at(threadNumber);
    unsigned lastFileIndex = UNMAPPED_BIN;
    FragmentBins lastFileBins;
    for (const unsigned binIndex : bins)
    {
        const unsigned fileIndex = binFiles_.at(binIndex);
        if (UNMAPPED_BIN != fileIndex)
        {
            if (lastFileIndex != fileIndex)
            {
                if (UNMAPPED_BIN != lastFileIndex)
                {
                    bufferBinIndexes(lastFileBins, buffers.at(lastFileIndex));
                    lastFileBins.clear();
                }

                FileBuffer &buffer = buffers[fileIndex];
                if (!bufferFragment(fragment, buffer))
                {
                    flushBuffer(buffer, binMetadataList, fileIndex);
                    ISAAC_VERIFY_MSG(bufferFragment(fragment, buffer), "Could not buffer into empty buffer" << fragment);
                }
                lastFileIndex = fileIndex;
            }

            lastFileBins.push_back(binIndex);
        }
    }
    bufferBinIndexes(lastFileBins, buffers.at(lastFileIndex));
}

void FragmentBinner::getFragmentStorageBins(const io::FragmentAccessor &fragment, FragmentBins &bins)
{
    if (!fragment.isAligned())
    {
        bins.push_back(0);
    }
    else
    {
        for (CigarPosition<const unsigned *> it(
            fragment.cigarBegin(), fragment.cigarEnd(), fragment.getFStrandReferencePosition(), fragment.isReverse(), fragment.readLength_);
            !it.end(); ++it)
        {
            if (Cigar::ALIGN == it.component().second)
            {
                const unsigned startBinIndex = binIndexMap_.getBinIndex(it.referencePos_);
                if (bins.empty() || bins.back() != startBinIndex)
                {
                    bins.push_back(startBinIndex);
                }
                // push the last position as well. This is important for duplicate detection of r-stranded alignments that end in the same bin but begin in different ones
                const unsigned endBinIndex = binIndexMap_.getBinIndex(it.referencePos_ + it.component().first - 1);
                if (bins.back() != endBinIndex)
                {
                    bins.push_back(endBinIndex);
                }
            }
        }
    }
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
}

void FragmentBinner::openBinFile(const BinMetadata &binMetadata, std::size_t file)
{
    ISAAC_THREAD_CERR << "openBin file: " << file << " for " << binMetadata << std::endl;
    // make sure file is empty first time we decide to put data in it.
    // boost::filesystem::remove for some stupid reason needs to allocate strings for this...
    if (common::deleteFile(binMetadata.getPath().c_str()) && ENOENT != errno)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to unlink " + binMetadata.getPath().string()));
    }
    files_[file].reopen(binMetadata.getPath().c_str(),
                        expectedBinSize_,
                        io::FileBufWithReopen::SequentialOnce);

    if (!files_[file].is_open())
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open bin file " + binMetadata.getPathString()));
    }
}

static std::size_t uniquePathCount(
    const alignment::BinMetadataList::const_iterator binsBegin,
    const alignment::BinMetadataList::const_iterator binsEnd)
{
    std::size_t ret = 0;
    std::unique_copy(
        boost::make_transform_iterator(binsBegin, boost::bind(&BinMetadata::getPath, _1)),
        boost::make_transform_iterator(binsEnd, boost::bind(&BinMetadata::getPath, _1)),
        boost::make_function_output_iterator([&ret](const boost::filesystem::path&){++ret;}));

    ISAAC_THREAD_CERR << "found " << ret << " unique bin paths in " << std::distance(binsBegin, binsEnd) << " bins";
    return ret;
}

void FragmentBinner::open(
    const BinMetadataList::iterator binsBegin,
    const BinMetadataList::iterator binsEnd)
{
    std::vector<io::FileBufWithReopen>(uniquePathCount(binsBegin, binsEnd), io::FileBufWithReopen(std::ios_base::out | std::ios_base::app | std::ios_base::binary)).swap(files_);
    binFiles_.resize(std::max_element(binsBegin, binsEnd, [](const BinMetadata& left, const BinMetadata& right){return left.getIndex() < right.getIndex();})->getIndex() + 1);
    
    ISAAC_TRACE_STAT("TemplateBuilder before Reopening output files");
    ISAAC_ASSERT_MSG(binsBegin != binsEnd, "Requested to open files for an empty bin range");
    ISAAC_ASSERT_MSG(std::size_t(std::distance(binsBegin, binsEnd)) <= binFiles_.size(),
                     "Requested more bins than anticipated during initialization:" << std::distance(binsBegin, binsEnd) << " max: " << binFiles_.size());
    ISAAC_THREAD_CERR << "Reopening output files for " << std::distance(binsBegin, binsEnd) << " bins" << std::endl;

    // if all needed files are open at the same time, there is no need to reopen anything.
    std::fill(binFiles_.begin(), binFiles_.begin() + binsBegin->getIndex(), UNMAPPED_BIN);
    std::fill(binFiles_.begin() + (binsEnd - 1)->getIndex() + 1, binFiles_.end(), UNMAPPED_BIN);


    alignment::BinMetadataList::iterator last = binsBegin;
    std::size_t file = 0;
    openBinFile(*binsBegin, file);
    for (alignment::BinMetadataList::iterator current = binsBegin; binsEnd != current; ++current)
    {
        // multiple BinMetadata may refer to the same storage file. Open each file only once
        if (last->getPath() != current->getPath())
        {
            ++file;
            openBinFile(*current, file);
        }
        binFiles_.at(current->getIndex()) = file;
//        ISAAC_THREAD_CERR << "mapped " << *current << " to file: " << file << std::endl;
        last = current;
    }

    for (FileBuffers &fileBuffers : threadFileBuffers_)
    {
        fileBuffers.clear();
        ISAAC_THREAD_CERR << "allocating " << files_.size() * sizeof(FileBuffer) << " bytes" << std::endl;
        fileBuffers.resize(files_.size());
    }

    ISAAC_THREAD_CERR << "Reopening output files done for " << std::distance(binsBegin, binsEnd) << " bins, reopened " << file << " files" << std::endl;
    ISAAC_TRACE_STAT("TemplateBuilder after Reopening output files");
}

void FragmentBinner::flush(BinMetadataList &binMetadataList)
{
    ISAAC_THREAD_CERR << "flushing " << files_.size() << " output buffers for " << threadFileBuffers_.size() << " threads "<< std::endl;
    for (FileBuffers &buffers : threadFileBuffers_)
    {
        unsigned fileIndex = 0;
        for (FileBuffer &buffer : buffers)
        {
            if (!buffer.empty())
            {
                flushBuffer(buffer, binMetadataList, fileIndex);
            }
            ++fileIndex;
        }
    }
    ISAAC_THREAD_CERR << "flushing " << files_.size() << " output buffers done for " << threadFileBuffers_.size() << " threads "<< std::endl;
}

void FragmentBinner::close() noexcept
{
    ISAAC_THREAD_CERR << "truncating " << files_.size() << " output files for " << std::endl;

    // make sure everything is written out for those that are open
    std::for_each(files_.begin(), files_.end(), boost::bind(&io::FileBufWithReopen::close, _1));

    std::fill(binFiles_.begin(), binFiles_.end(), UNMAPPED_BIN);

    ISAAC_THREAD_CERR << "truncating done for " << files_.size() << " output files" << std::endl;
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
