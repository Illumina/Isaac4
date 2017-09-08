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
 ** \file Contig.
 **
 ** \brief See Contig.hh
 **
 ** \author Come Raczy
 **/

#include <algorithm>
#include <numeric>
#include <boost/bind.hpp>

#include "reference/Contig.hh"

namespace isaac
{
namespace reference
{

static std::size_t roundToPadding(const std::size_t size, const std::size_t padding)
{
    return (((size + padding - 1) / padding) * padding);
}

static std::size_t genomeSize(
    const isaac::reference::SortedReferenceMetadata::Contigs &contigList,
    const std::size_t padding,
    const std::size_t spacing)
{
    return std::accumulate(
        contigList.begin(), contigList.end(),
        size_t(0), [padding, spacing](const std::size_t sum, const SortedReferenceMetadata::Contig &contig)
                           {   return sum + roundToPadding(contig.totalBases_ + spacing, padding);}) + spacing;
}

/**
 * \brief construct a reference memory block with contigs placed so that there is at least
 *        spacing number of bytes between them and spacing number of bytes after the last contig
 *
 * @param spacing   The minimum number of bytes allocated in front of each contig and after the last contig that are not part of any contig.
 *                  This allows for fast and simple mismatch counting that does not need to consider the edge cases
 * @param padding   No two contigs will produce collision when any of their absolute positions is divided by padding. This way
 *                  the memory offsets can be unambiguously mapped back to contig index with a simple division and a lookup.
 */
template <typename AllocatorT>
BasicContigList<AllocatorT>::BasicContigList(
    const isaac::reference::SortedReferenceMetadata::Contigs &contigMetadataList,
    const std::size_t spacing):
    contigIdFromScaledOffset_(TRANSLATION_TABLE_SIZE, INVALID_CONTIG_ID),
    referenceSequence_(genomeSize(contigMetadataList, CONTIG_LENGTH_MIN, spacing))
{
    this->reserve(contigMetadataList.size() + 1);
    ReferenceSequenceIterator b = referenceSequence_.begin();
    std::transform(
        contigMetadataList.begin(), contigMetadataList.end(), std::back_inserter(base()),
        [this, &b, spacing](const isaac::reference::SortedReferenceMetadata::Contig &contigMetadata)
        {
            Contig ret(contigMetadata.index_, contigMetadata.name_, contigMetadata.decoy_, b + spacing, b + spacing + contigMetadata.totalBases_);
            b += roundToPadding(contigMetadata.totalBases_ + spacing, CONTIG_LENGTH_MIN);
            ISAAC_ASSERT_MSG(referenceSequence_.end() >= b, "Overrun");
            return ret;
        });
    // the fake one which is there just to provide spacing for the last real one
    this->push_back(Contig(-1, "", false, b + spacing, b + spacing));

    for (ContigId contigId = 0; contigId < size() - 1; ++contigId)
    {
        for (std::size_t offset = contigBeginOffset(contigId); offset < contigBeginOffset(contigId + 1); offset += CONTIG_LENGTH_MIN)
        {
            contigIdFromScaledOffset_[offset / CONTIG_LENGTH_MIN] = contigId;
        }
    }
    // remove fake one to avoid messing up end();
    this->pop_back();

    const std::size_t decoys =
        std::count_if(begin(), end(), [](const Contig &contig){return contig.isDecoy();});
    ISAAC_THREAD_CERR << "Generated " << size() << " contigs of which " << decoys << " are decoys" << std::endl;

}

template class BasicContigList<common::NumaAllocator<char, 0> >;

} // namespace reference
} // namespace isaac
