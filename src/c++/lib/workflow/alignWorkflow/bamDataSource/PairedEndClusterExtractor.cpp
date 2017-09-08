/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file PairedEndClusterExtractor.cpp
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/integer/static_min_max.hpp>

#include "common/Debug.hh"
#include "workflow/alignWorkflow/bamDataSource/PairedEndClusterExtractor.hh"
#include "oligo/Nucleotides.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace bamDataSource
{

void TempFileClusterExtractor::open(const boost::filesystem::path &tempFilePath, std::streamsize expectedFileSize)
{
    if (tempFilePath_ != tempFilePath)
    {
        recordIndex_.clear();
        unpairedReadsFile_.reopen(tempFilePath.c_str(), io::FileBufWithReopen::SequentialOnce);
        std::istream is(&unpairedReadsFile_);
        if (!is.eof())
        {
            if (!is)
            {
                BOOST_THROW_EXCEPTION(isaac::common::IoException(
                    errno, (boost::format("Unable to open file for reading %s") % tempFilePath.string()).str()));
            }

            resize(expectedFileSize + 1);
            if (!is.read(&front(), size()) && !is.eof())
            {
                BOOST_THROW_EXCEPTION(isaac::common::IoException(
                    errno, (boost::format("Unable to read %d bytes from %s") % size() % tempFilePath.string()).str()));
            }

            ISAAC_ASSERT_MSG(is.gcount() == expectedFileSize, "Read mismatching number of bytes from the file: " << expectedFileSize <<
                " from " << tempFilePath);

            resize(expectedFileSize);

            for (std::vector<char>::const_iterator it = begin(); end() != it; it = getNextRecord(it))
            {
                recordIndex_.push_back(it);
            }

            std::sort(recordIndex_.begin(), recordIndex_.end(), compareNameAndRead);
            ISAAC_THREAD_CERR << "TempFileClusterExtractor::open: " << recordIndex_.size() << " " << tempFilePath << std::endl;
        }
        else
        {
            ISAAC_THREAD_CERR << "TempFileClusterExtractor::open: empty " << tempFilePath << std::endl;
        }
    }
    firstUnextracted_ = recordIndex_.begin();
    tempFilePath_ = tempFilePath.c_str();
}

template <typename IteratorT>
void UnpairedReadsCache::storeUnpaired(
    IteratorT unpairedBegin,
    IteratorT unpairedEnd,
    const flowcell::ReadMetadataList &readMetadataList)
{
#ifdef ISAAC_DEV_STATS_ENABLED
//        ISAAC_ASSERT_MSG(false, "Unpaired resolution is not supported with ISAAC_DEV_STATS_ENABLED");
    return;
#endif// ISAAC_DEV_STATS_ENABLED
    common::StaticVector<char, 10240> byteBuff;
    BOOST_FOREACH(const IndexRecord &idx, std::make_pair(unpairedBegin, unpairedEnd))
    {
        const flowcell::ReadMetadata &readMetadata = readMetadataList.at(!idx.getBlock().isReadOne());

        const bam::BamBlockHeader &block = idx.getBlock();

        byteBuff.clear();
        extractReadName(block, maxReadNameLength_, std::back_inserter(byteBuff));
        ISAAC_ASSERT_MSG(byteBuff.size() == maxReadNameLength_, "Invalid number of name bytes extracted");
        byteBuff.push_back(0);// 0 terminator is needed for name comparison during extraction
        const unsigned nameCrc = getNameCrc<7>(crcWidth_, byteBuff.begin(), maxReadNameLength_);
        std::ostream os(tempFiles_[nameCrc].get());

        const TempFileClusterExtractor::FlagsType flags =
            (block.isReadOne() ? TempFileClusterExtractor::READ_ONE_FLAG : 0) |
                (block.isPf() ? TempFileClusterExtractor::PASS_FILTER_FLAG : 0) ;
        const unsigned recordLength = sizeof(recordLength) + sizeof(flags) + byteBuff.size() + readMetadata.getLength();
        if (!os.write(reinterpret_cast<const char*>(&recordLength), sizeof(unsigned)))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % sizeof(unsigned) % common::pathStringToStdString(tempFilePaths_[nameCrc])).str()));
        }
        tempFileSizes_[nameCrc] += sizeof(unsigned);

        if (!os.write(reinterpret_cast<const char*>(&flags), sizeof(flags)))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % sizeof(flags) % common::pathStringToStdString(tempFilePaths_[nameCrc])).str()));
        }
        tempFileSizes_[nameCrc] += sizeof(flags);

        if (!os.write(byteBuff.begin(), byteBuff.size()))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % byteBuff.size() % common::pathStringToStdString(tempFilePaths_[nameCrc])).str()));
        }
        tempFileSizes_[nameCrc] += byteBuff.size();

        byteBuff.resize(readMetadata.getLength());
        bam::extractBcl(idx.getBlock(), byteBuff.begin(), readMetadata);

        if (!os.write(&byteBuff.front(), byteBuff.size()))
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(
                errno, (boost::format("Failed to write: %d bytes into %s") % byteBuff.size() % common::pathStringToStdString(tempFilePaths_[nameCrc])).str()));
        }
        tempFileSizes_[nameCrc] += byteBuff.size();
    }
}

void PairedEndClusterExtractor::reset()
{
    firstUnextracted_ = end();
}

inline bool inRange(const IndexRecord &idx, const char* rangeStart, const char* rangeEnd)
{
    return reinterpret_cast<const char*>(idx.bamRecordPointer_) >= rangeStart &&
        reinterpret_cast<const char*>(idx.bamRecordPointer_) < rangeEnd;
}

/**
 * \brief stores bam records in [rangeStart,rangeEnd) into temporary files and frees memory associated with them
 * \precondition All records in the extractor are assumed to be unpaired
 */
void PairedEndClusterExtractor::removeOld(
    const char* rangeStart,
    const char* rangeEnd,
    const flowcell::ReadMetadataList &readMetadataList)
{
    BaseT::iterator oldBegin =
        std::partition(begin(), end(), !boost::bind(&inRange, _1, rangeStart, rangeEnd));

//    ISAAC_THREAD_CERR << "remove old " << std::size_t(rangeStart) << "-" << std::size_t(rangeEnd) << " oldBegin!=begin():" << (oldBegin != begin()) << std::endl;
    storeUnpaired(oldBegin, end(), readMetadataList);

    erase(oldBegin, end());
    reset();
}

void PairedEndClusterExtractor::sort()
{
    std::sort(begin(), end());
    firstUnextracted_ = begin();
}

template
unsigned PairedEndClusterExtractor::extractUnpaired<std::vector<char>::iterator, std::vector<bool>::iterator >(
    const unsigned r1Length,
    const unsigned r2Length,
    const unsigned nameLengthMax,
    unsigned clusterCount, std::vector<char>::iterator &clusterIt, std::vector<bool>::iterator &pfIt );

template
unsigned PairedEndClusterExtractor::extractPairedReads<std::vector<char>::iterator, std::vector<bool>::iterator >(
    const unsigned nameLengthMax, unsigned clusterCount, std::vector<char>::iterator &clusterIt, std::vector<bool>::iterator &pfIt,
    const flowcell::ReadMetadataList &readMetadataList);

} // namespace bamDataSource
} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
