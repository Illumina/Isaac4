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
 ** \file SemialignedEndsClipper.cpp
 **
 ** \brief See SemialignedEndsClipper.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/Mismatch.hh"
#include "common/Debug.hh"

#include "SemialignedEndsClipper.hh"

namespace isaac
{
namespace build
{

using alignment::Cigar;

/**
 * \brief clips mismatches on the left if this does not move index.pos_ to binEndPos or beyond
 */
bool SemialignedEndsClipper::clipLeftSide(
    const reference::ContigList &contigList,
    const reference::ReferencePosition binEndPos,
    const unsigned char *sequenceBegin,
    PackedFragmentBuffer::Index &index,
    unsigned short &newEditDistance)
{
    const uint32_t *oldCigarBegin = index.cigarBegin_;
    std::pair<unsigned, alignment::Cigar::OpCode> operation = alignment::Cigar::decode(*oldCigarBegin);
    unsigned softClippedBeginBases = 0;
    if (alignment::Cigar::SOFT_CLIP == operation.second)
    {
        ++oldCigarBegin;
        softClippedBeginBases = operation.first;
        sequenceBegin += operation.first;
        operation = alignment::Cigar::decode(*oldCigarBegin);
    }

    if (alignment::Cigar::ALIGN == operation.second)
    {
        unsigned mappedBeginBases = operation.first;
        const unsigned char * sequenceEnd = sequenceBegin + mappedBeginBases;

        const reference::Contig &contig = contigList.at(index.pos_.getContigId());
        reference::Contig::const_iterator referenceBegin = contig.begin() + index.pos_.getPosition();

        std::pair<unsigned, unsigned> clipped = alignment::clipMismatches<CONSECUTIVE_MATCHES_MIN>(sequenceBegin, sequenceEnd,
                                                                              referenceBegin, contig.end(),
                                                                              &oligo::getReferenceBaseFromBcl);

        if (clipped.first && index.pos_ + clipped.first < binEndPos)
        {
            softClippedBeginBases += clipped.first;
            mappedBeginBases -= clipped.first;
            index.pos_ += clipped.first;
            // don't update fStrandPosition_ here. index has it and it will be properly synced by gap realigner updatePairDetails
//            fragment.fStrandPosition_ += clipped.first;
            newEditDistance -= clipped.second;

            size_t before = cigarBuffer_.size();
            cigarBuffer_.addOperation(softClippedBeginBases, Cigar::SOFT_CLIP);
            cigarBuffer_.addOperation(mappedBeginBases, Cigar::ALIGN);
            cigarBuffer_.addOperations(oldCigarBegin + 1, index.cigarEnd_);
            index.cigarBegin_ = &cigarBuffer_.at(before);
            index.cigarEnd_ = &cigarBuffer_.back() + 1;
            return true;
        }
    }
    return false;
}

bool SemialignedEndsClipper::clipRightSide(
    const reference::ContigList &contigList,
    const unsigned char *sequenceEnd,
    PackedFragmentBuffer::Index &index,
    reference::ReferencePosition &newRStrandPosition,
    unsigned short &newEditDistance)
{
    std::reverse_iterator<const unsigned char *> sequenceRBegin(sequenceEnd);

    const uint32_t *oldCigarEnd = index.cigarEnd_;
    std::pair<unsigned, alignment::Cigar::OpCode> operation = alignment::Cigar::decode(*(oldCigarEnd - 1));
    unsigned softClippedEndBases = 0;
    if (Cigar::SOFT_CLIP == operation.second)
    {
        --oldCigarEnd;
        softClippedEndBases = operation.first;
        sequenceRBegin += operation.first;
        operation = alignment::Cigar::decode(*(oldCigarEnd - 1));
    }

    if (alignment::Cigar::ALIGN == operation.second)
    {
        unsigned mappedEndBases = operation.first;
        std::reverse_iterator<const unsigned char *> sequenceREnd = sequenceRBegin + mappedEndBases;

        const reference::Contig &reference = contigList.at(index.pos_.getContigId());
//        ISAAC_ASSERT_MSG(index.pos_ == fragment.fStrandPosition_, "Broken index :" << index << "\n" << fragment);
        std::reverse_iterator<reference::Contig::const_iterator> referenceRBegin(reference.begin() + newRStrandPosition.getPosition() + 1);
        std::reverse_iterator<reference::Contig::const_iterator> referenceREnd(reference.begin());

        std::pair<unsigned, unsigned> clipped = alignment::clipMismatches<CONSECUTIVE_MATCHES_MIN>(sequenceRBegin, sequenceREnd,
                                                                              referenceRBegin, referenceREnd,
                                                                              &oligo::getReferenceBaseFromBcl);

        if (clipped.first)
        {
            softClippedEndBases += clipped.first;
            mappedEndBases -= clipped.first;
            newRStrandPosition -= clipped.first;
            newEditDistance -= clipped.second;

            size_t before = cigarBuffer_.size();
            cigarBuffer_.addOperations(index.cigarBegin_, oldCigarEnd - 1);
            cigarBuffer_.addOperation(mappedEndBases, Cigar::ALIGN);
            cigarBuffer_.addOperation(softClippedEndBases, Cigar::SOFT_CLIP);
            index.cigarBegin_ = &cigarBuffer_.at(before);
            index.cigarEnd_ = &cigarBuffer_.back() + 1;
            return true;
        }
    }
    return false;
}

void SemialignedEndsClipper::clip(
    const reference::ContigList &contigs,
    const reference::ReferencePosition binEndPos,
    const io::FragmentAccessor &fragment,
    PackedFragmentBuffer::Index &index,
    reference::ReferencePosition &newRStrandPosition,
    unsigned short &newEditDistance)
{
    ISAAC_ASSERT_MSG(fragment.isAligned(), "Unexpected unaligned fragment from gap realigner");

    const bool leftClipped = clipLeftSide(contigs, binEndPos, fragment.basesBegin(), index, newEditDistance);
    const bool rightClipped = clipRightSide(contigs, fragment.basesEnd(), index, newRStrandPosition, newEditDistance);
    if (leftClipped || rightClipped)
    {
        ISAAC_THREAD_CERR_DEV_TRACE(" SemialignedEndsClipper::clip: " << fragment);
    }
}

} // namespace build
} // namespace isaac
