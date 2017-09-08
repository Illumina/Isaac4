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
 ** \file MergeReferenceWorkflow.cpp
 **
 ** \brief see MergeReferenceWorkflow.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "workflow/MergeReferencesWorkflow.hh"

namespace isaac
{
namespace workflow
{

MergeReferencesWorkflow::MergeReferencesWorkflow(
    const std::vector<bfs::path> &filesToMerge,
    const bfs::path &outputFilePath,
    const bool makeAbsolutePaths)
    : filesToMerge_(filesToMerge),
      outputFilePath_(outputFilePath),
      makeAbsolutePaths_(makeAbsolutePaths)
{
}

void MergeReferencesWorkflow::run()
{
    reference::SortedReferenceMetadata result;
    BOOST_FOREACH(const boost::filesystem::path &path, filesToMerge_)
    {
        reference::SortedReferenceMetadata referenceToMerge = reference::loadReferenceMetadataFromXml(path, makeAbsolutePaths_);
        result.merge(referenceToMerge);
    }
    const reference::SortedReferenceMetadata::Contigs karyotypeOrderedContigs = result.getContigs();
    const reference::SortedReferenceMetadata::Contigs::const_iterator collision = std::adjacent_find(
        karyotypeOrderedContigs.begin(), karyotypeOrderedContigs.end(),
        boost::bind(&reference::SortedReferenceMetadata::Contig::index_, _1) ==
            boost::bind(&reference::SortedReferenceMetadata::Contig::index_, _2));

    if (karyotypeOrderedContigs.end() != collision)
    {
        const boost::format message = boost::format("\n   *** Karyotype index collision detected in %s and %s ***\n") % *collision % *(collision + 1);
        BOOST_THROW_EXCEPTION(common::PostConditionException(message.str()));
    }

    reference::saveSortedReferenceXml(outputFilePath_, result);
}

} // namespace workflow
} // namespace isaac
