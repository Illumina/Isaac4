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
 ** \file SortedReferenceXml.hh
 **
 ** sorted-reference.xml io helper.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_SORTED_REFERENCE_FASTA_HH
#define iSAAC_REFERENCE_SORTED_REFERENCE_FASTA_HH

#include <boost/filesystem.hpp>

#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{
SortedReferenceMetadata loadReferenceMetadataFromFasta(
    const boost::filesystem::path &xmlPath,
    common::ThreadVector &threads);

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_SORTED_REFERENCE_FASTA_HH
