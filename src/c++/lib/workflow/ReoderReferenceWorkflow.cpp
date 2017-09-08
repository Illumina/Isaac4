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
 ** \file ReoderReferenceWorkflow.cpp
 **
 ** \brief see ReoderReferenceWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/algorithm/string/join.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/ContigLoader.hh"
#include "workflow/ReorderReferenceWorkflow.hh"

namespace isaac
{
namespace workflow
{

ReorderReferenceWorkflow::ReorderReferenceWorkflow(
    const bfs::path &sortedReferenceMetadata,
    const bfs::path &newXmlPath,
    const bfs::path &newDataFileDirectory,
    const std::vector<std::string> &newOrder
    )
    : sortedReferenceMetadata_(sortedReferenceMetadata),
      newXmlPath_(newXmlPath),
      newDataFileDirectory_(newDataFileDirectory),
      xml_(reference::loadReferenceMetadataFromXml(sortedReferenceMetadata_))
{
    const reference::SortedReferenceMetadata::Contigs &contigs = xml_.getContigs();
    if (newOrder.empty())
    {
        originalIndexes_.clear();
        std::transform(contigs.begin(), contigs.end(), std::back_inserter(originalIndexes_),
                       boost::bind(&reference::SortedReferenceMetadata::Contig::index_, _1));

        ISAAC_THREAD_CERR << "Preserving the existing order of contigs" << std::endl;
        return;
    }

    typedef std::pair<std::string, unsigned> NameIndex;

    std::vector<NameIndex > newOrderIndexed;
    newOrderIndexed.reserve(newOrder.size());
    for (const std::string &name : newOrder)
    {
        newOrderIndexed.push_back(std::make_pair(name, newOrderIndexed.size()));
    }
    std::sort(newOrderIndexed.begin(), newOrderIndexed.end(), [](const NameIndex &left, const NameIndex &right){return left.first < right.first;});

    std::vector<NameIndex > oldOrderIndexed;
    oldOrderIndexed.reserve(contigs.size());
    for (const reference::SortedReferenceMetadata::Contig &xmlContig : contigs)
    {
        oldOrderIndexed.push_back(std::make_pair(xmlContig.name_, oldOrderIndexed.size()));
    }
    std::sort(oldOrderIndexed.begin(), oldOrderIndexed.end(), [](const NameIndex &left, const NameIndex &right){return left.first < right.first;});


    std::vector<unsigned> to(contigs.size());
    std::iota(to.begin(), to.end(), 0);
    std::vector<unsigned> from(to);

    // move elements to their new indexes while trying not to touch
    // the ones that are not mentioned in newOrder
    for (auto o = oldOrderIndexed.begin(), n = newOrderIndexed.begin();
            oldOrderIndexed.end() != o && newOrderIndexed.end() != n;
            )
    {
        if (o->first < n->first)
        {
            ++o;
        }
        else if (n->first < o->first)
        {
            ++n;
        }
        else
        {
            unsigned oi = o->second;
            unsigned ni = n->second;

            to[from[ni]] = to[oi];
            from[to[oi]] = from[ni];

            to[oi] = ni;
            from[ni] = oi;
            ++o;
            ++n;
        }
    }

    originalIndexes_.resize(from.size());

    for (std::size_t i = 0; i < originalIndexes_.size(); ++i)
    {
        originalIndexes_[i] = from[i];
    }
}

void ReorderReferenceWorkflow::run()
{
    std::ofstream xmlOs(newXmlPath_.c_str());
    if (!xmlOs)
    {
        BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open output file: " + newXmlPath_.string()));
    }

    const boost::filesystem::path targetPath = newDataFileDirectory_ / xml_.getContigs().front().filePath_.filename();
    std::ofstream ofs(targetPath.c_str(), std::ios_base::binary);
    std::size_t newIndex = 0;
    for(unsigned idx: originalIndexes_)
    {
        reference::SortedReferenceMetadata::Contig &xmlContig = xml_.getContigs().at(idx);
        std::ifstream ifs(xmlContig.filePath_.c_str(), std::ios_base::binary);
        if (!ifs)
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to open input file: " + xmlContig.filePath_.string()));
        }
        ifs.seekg(xmlContig.offset_);
        if (!ifs)
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to seek to position " + std::to_string(xmlContig.offset_) +
                " in the input file: " + xmlContig.filePath_.string()));
        }
        ofs << ">" << xmlContig.name_ << std::endl;
        if (!ifs)
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to write contig name " + xmlContig.name_ +
                " into the output file: " + xmlContig.filePath_.string()));
        }
        xmlContig.offset_ = ofs.tellp();
        std::copy_n(std::istreambuf_iterator<char>(ifs), xmlContig.size_, std::ostreambuf_iterator<char>(ofs));
        if (!ifs)
        {
            BOOST_THROW_EXCEPTION(isaac::common::IoException(errno, "Failed to copy" + std::to_string(xmlContig.size_) +
                " bytes from input file: " + xmlContig.filePath_.string() +
                " into the output file: " + targetPath.string()));
        }
        xmlContig.filePath_ = targetPath;
        xmlContig.index_ = newIndex++;
        ISAAC_THREAD_CERR << "Copied " << xmlContig << std::endl;
    }
    std::sort(xml_.getContigs().begin(), xml_.getContigs().end(),
            [](const reference::SortedReferenceMetadata::Contig &left,
                    const reference::SortedReferenceMetadata::Contig &right){return left.index_ < right.index_;});

    saveSortedReferenceXml(xmlOs, xml_);
}

} // namespace workflow
} // namespace isaac
