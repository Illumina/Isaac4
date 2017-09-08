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
 ** \file SortedReferenceXml.cpp
 **
 ** SortedReference.xml helper.
 **
 ** \author Roman Petrovski
 **/

#include <mutex>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/MD5Sum.hh"
#include "common/SystemCompatibility.hh"
#include "common/Threads.hpp"
#include "reference/SortedReferenceFasta.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

static reference::SortedReferenceMetadata::Contig parseContig(
    std::vector<char> &fileContents,
    std::vector<char>::iterator begin,
    std::vector<char>::iterator end)
{
    const std::string contigName(
        begin + 1, std::find_if(begin + 1, end, [](const char c){return ' ' == c || '\n' == c || '\r' == c || '\t' == c;}));

    const std::vector<char>::iterator firstNewline = std::find_if(begin, end, [](char c){return c == '\r' || c == '\n';});
    const std::vector<char>::iterator firstBase = std::find_if(firstNewline, end, [](char c){return c != '\r' && c != '\n';});

    const std::size_t acgtCount = std::count_if(firstBase, end, [](char c){return oligo::getValue(c) != oligo::INVALID_OLIGO;});
    const std::size_t totalBasesCount = std::count_if(firstBase, end, [](char c){return c != '\r' && c != '\n';});

    isaac::common::MD5Sum md5Sum;

    const std::size_t byteSize = std::distance(firstBase, end);
    end = std::remove_if(firstBase, end, [](const char c){return ' ' == c || '\n' == c || '\r' == c || '\t' == c;});
    auto range = boost::make_iterator_range(firstBase, end);
    boost::to_upper(range);

    md5Sum.update(&*firstBase, std::distance(firstBase, end));

    reference::SortedReferenceMetadata::Contig ret(
        -1,
        contigName,
        false,
        "",
        std::distance(fileContents.begin(), firstBase),
        byteSize,
        -1,
        totalBasesCount,
        acgtCount,
        "",
        "",
        isaac::common::MD5Sum::toHexString( md5Sum.getDigest().data, 16 ));

    return ret;
}

static void parseFastaThread(
    std::vector<char> &fileContents,
    std::mutex &m,
    std::size_t &offset,
    std::size_t &index,
    const boost::filesystem::path &fastaPath,
    SortedReferenceMetadata &sortedReferenceMetadata
    )
{
    std::lock_guard<std::mutex> lock(m);

    while (offset != fileContents.size())
    {
        ISAAC_ASSERT_MSG('>' == fileContents.at(offset), "Offset " << offset << " does not point at '>'");
        const std::vector<char>::iterator nextHeader = std::find(fileContents.begin() + offset + 1, fileContents.end(), '>');
        const std::size_t beginOffset = offset;
        offset = std::distance(fileContents.begin(), nextHeader);
        const std::size_t endOffset = offset;
        const std::size_t ourIndex = index;
        ++index;

        reference::SortedReferenceMetadata::Contig contig;
        {
            common::unlock_guard<std::mutex > unlock(m);
            contig = parseContig(fileContents, fileContents.begin() + beginOffset, fileContents.begin() + endOffset);
        }
        contig.index_ = ourIndex;
        contig.filePath_ = fastaPath;
        ISAAC_THREAD_CERR << "found " << contig << std::endl;
        sortedReferenceMetadata.putContig(contig);
    }
}

SortedReferenceMetadata loadReferenceMetadataFromFasta(
    const boost::filesystem::path &fastaPath,
    common::ThreadVector &threads)
{
    std::vector<char> fileContents(common::getFileSize(fastaPath.c_str()));
    std::ifstream is(fastaPath.c_str());
    if (!is)
    {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open sorted reference file " + fastaPath.string()));
    }

    //load entire file in memory.
    is.read(&fileContents.front(), fileContents.size());
    if (!is && !is.eof())
    {
        BOOST_THROW_EXCEPTION(common::IoException(
            errno, (boost::format("Failed read %d bytes from %s. Read %d") % fileContents.size() % fastaPath.string() % is.tellg()).str()));
    }

    if (fileContents.empty() || '>' != fileContents.front())
    {
        BOOST_THROW_EXCEPTION(common::IoException(
            0, (boost::format("Fasta file is expected to have '>' as first character: %s") % fastaPath.string()).str()));
    }

    SortedReferenceMetadata ret;
    std::mutex m;
    std::size_t offset = 0;
    std::size_t index = 0;
    threads.execute([&fileContents, &m, &offset, &index, &fastaPath, &ret](unsigned threadNumber, unsigned threadsTotal)
                    {parseFastaThread(fileContents, m, offset, index, fastaPath, ret);});

    reference::SortedReferenceMetadata::Contigs &contigs = ret.getContigs();
    std::sort(contigs.begin(), contigs.end(),
              [](const reference::SortedReferenceMetadata::Contig &left,
                  const reference::SortedReferenceMetadata::Contig &right){return left.index_ < right.index_;});
    std::size_t genomicOffset = 0;
    for (reference::SortedReferenceMetadata::Contig &contig : contigs)
    {
        contig.genomicPosition_ = genomicOffset;
        genomicOffset += contig.totalBases_;
    }

    ret.makeAbsolutePaths(fastaPath.parent_path());

    return ret;
}

} // namespace reference
} // namespace isaac


