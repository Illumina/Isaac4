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
 ** \file PeAdapterTrimmer.cpp
 **
 ** \author Roman Petrovski
 **/

#include "alignment/templateBuilder/PeAdapterTrimmer.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

PeAdapterTrimmer::PeAdapterTrimmer(
    const bool collectMismatchCycles,
    const bool trimPEAdapters,
    const AlignmentCfg &alignmentCfg)
    : trimPEAdapters_(trimPEAdapters)
    , collectMismatchCycles_(collectMismatchCycles)
    , alignmentCfg_(alignmentCfg)
{
}

FragmentMetadata PeAdapterTrimmer::trimForwardPEAdapter(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const FragmentMetadata & forwardFragment,
    const reference::ReferencePosition &adapterPosition,
    Cigar &cigarBuffer) const
{
    // trim forward fragment
    const unsigned cigarOffset = cigarBuffer.size();
    // ensure reallocation does not occur during the copying.
//    cigarBuffer_.reserve(cigarBuffer_.size() + std::distance(forwardFragment.cigarBegin(), forwardFragment.cigarEnd()));
    cigarBuffer.addOperations(forwardFragment.cigarBegin(), forwardFragment.cigarEnd());
    reference::ReferencePosition revPos = forwardFragment.getRStrandReferencePosition();
    uint64_t clip = 0;
    uint64_t align = 0;
    while (cigarOffset != cigarBuffer.size())
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer.back());
        cigarBuffer.pop_back();
        if(Cigar::SOFT_CLIP == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::INSERT == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::DELETE == cigar.second)
        {
            const int64_t rightOverhang = revPos - adapterPosition;
            if (rightOverhang < cigar.first)
            {
                break;
            }
            revPos -= cigar.first;
        }
        else if (Cigar::ALIGN == cigar.second)
        {
            const int64_t rightOverhang = revPos - adapterPosition;
            if (rightOverhang < cigar.first)
            {
                clip += rightOverhang;
                align = cigar.first - rightOverhang;
                break;
            }
            clip += cigar.first;
            revPos -= cigar.first;
        }
        else
        {
            ISAAC_VERIFY_MSG(false, "Unexpected cigar operation:" << cigar << " " << forwardFragment);
        }
    }
    if (align)
    {
        cigarBuffer.addOperation(align, Cigar::ALIGN);
    }
    if (clip)
    {
        cigarBuffer.addOperation(clip, Cigar::SOFT_CLIP);
    }
    FragmentMetadata ret = forwardFragment;
    ret.resetAlignment();
    ret.rightClipped() = std::max<unsigned short>(clip, ret.rightClipped());
    ret.incrementAdapterClip(clip);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " before updateAlignment:" << ret);

    ret.updateAlignment(
        collectMismatchCycles_, alignmentCfg_, readMetadata, contigList,
        forwardFragment.reverse,
        forwardFragment.contigId, forwardFragment.position,
        cigarBuffer, cigarOffset);
    return ret;
}

bool PeAdapterTrimmer::trimForwardPEAdapter(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    uint32_t cycles,
    FragmentMetadata & forwardFragment,
    Cigar &cigarBuffer) const
{
    if (cycles <= forwardFragment.rightClipped())
    {
        return false;
    }

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " trim fw:" << forwardFragment << " cycles:" << cycles);

    const unsigned cigarOffset = cigarBuffer.size();
    cigarBuffer.addOperations(forwardFragment.cigarBegin(), forwardFragment.cigarEnd());
    uint64_t clip = 0;
    uint64_t align = 0;
    while (cigarOffset != cigarBuffer.size())
    {
        const Cigar::Component cigar = Cigar::decode(cigarBuffer.back());
        cigarBuffer.pop_back();
        if(Cigar::SOFT_CLIP == cigar.second)
        {
            clip += cigar.first;
            cycles -= std::min(cycles, cigar.first);
        }
        else if(Cigar::HARD_CLIP == cigar.second)
        {
            // Currently hard clips are only present in inversions. Too complicated to adjust, just bail out
            return false;
        }
        else if (Cigar::INSERT == cigar.second)
        {
            clip += cigar.first;
            cycles -= std::min(cycles, cigar.first);
        }
        else if (Cigar::DELETE == cigar.second || Cigar::BACK == cigar.second)
        {
            //these don't contain cycles, just remove them from cigar.
        }
        else if (Cigar::FLIP == cigar.second || Cigar::CONTIG == cigar.second)
        {
            //these are too complicated to adjust, just bail.
            return false;
        }
        else if (Cigar::ALIGN == cigar.second)
        {
            if (cycles < cigar.first)
            {
                clip += cycles;
                align = cigar.first - cycles;
                cycles = 0;
                break;
            }
            clip += cigar.first;
            cycles -= cigar.first;
        }
        else
        {
            ISAAC_VERIFY_MSG(false, "Unexpected cigar operation:" << cigar << " " << forwardFragment);
        }
    }

    // save the position as it changes during reset
    const int64_t pos = forwardFragment.position;
    forwardFragment.resetAlignment();
    if (!align && cigarBuffer.size() == cigarOffset)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " fultrim:" << forwardFragment);
    }
    else
    {
        ISAAC_ASSERT_MSG(!cycles, "Some cycles left unclipped:" << cycles << " " << forwardFragment);

        if (align)
        {
            cigarBuffer.addOperation(align, Cigar::ALIGN);
        }
        if (clip)
        {
            cigarBuffer.addOperation(clip, Cigar::SOFT_CLIP);
        }
        forwardFragment.rightClipped() = std::max<unsigned short>(clip, forwardFragment.rightClipped());
        forwardFragment.incrementAdapterClip(clip);
        forwardFragment.updateAlignment(
            collectMismatchCycles_, alignmentCfg_, readMetadata, contigList,
            forwardFragment.reverse,
            forwardFragment.contigId, pos,
            cigarBuffer, cigarOffset);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " trimmed:" << forwardFragment);
    }
    return true;
}

FragmentMetadata PeAdapterTrimmer::trimReversePEAdapter(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const FragmentMetadata & reverseFragment,
    const reference::ReferencePosition &adapterPosition,
    Cigar &cigarBuffer) const
{
    const unsigned cigarOffset = cigarBuffer.size();
    uint64_t clip = 0;
    uint64_t align = 0;
    reference::ReferencePosition fwPos = reverseFragment.getFStrandReferencePosition();

    Cigar::const_iterator current = reverseFragment.cigarBegin();
    for(;
        reverseFragment.cigarEnd() != current;
        ++current)
    {
        const Cigar::Component cigar = Cigar::decode(*current);
        if(Cigar::SOFT_CLIP == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::INSERT == cigar.second)
        {
            clip += cigar.first;
        }
        else if (Cigar::DELETE == cigar.second)
        {
            const int64_t leftOverhang = adapterPosition - fwPos;
            if (leftOverhang < cigar.first)
            {
                break;
            }
            fwPos += cigar.first;
        }
        else if (Cigar::ALIGN == cigar.second)
        {
            const int64_t leftOverhang = adapterPosition - fwPos;
            if (leftOverhang < cigar.first)
            {
                clip += leftOverhang;
                align = cigar.first - leftOverhang;
                fwPos += leftOverhang;
                break;
            }
            clip += cigar.first;
            fwPos += cigar.first;
        }
        else
        {
            ISAAC_VERIFY_MSG(false, "Unexpected cigar operation:" << cigar);
        }
    }

    ISAAC_ASSERT_MSG(reverseFragment.cigarEnd() != current, "Unexpectedly clipped everything for " << reverseFragment << " with " << adapterPosition);
    // ensure reallocation does not occur during the copying.
//    cigarBuffer_.reserve(cigarBuffer_.size() + reverseFragment.cigarLength - currentOffset + 2);
    if (clip)
    {
        cigarBuffer.addOperation(clip, Cigar::SOFT_CLIP);
    }
    if (align)
    {
        cigarBuffer.addOperation(align, Cigar::ALIGN);
    }
    cigarBuffer.addOperations(current + 1, reverseFragment.cigarEnd());

    FragmentMetadata ret = reverseFragment;
    ret.resetAlignment();
    ret.updateAlignment(
        collectMismatchCycles_, alignmentCfg_, readMetadata, contigList,
        reverseFragment.reverse,
        reverseFragment.contigId, fwPos.getPosition(),
        cigarBuffer, cigarOffset);
    ret.leftClipped() = std::max<unsigned short>(clip, ret.leftClipped());
    ret.incrementAdapterClip(clip);
    return ret;
}

bool PeAdapterTrimmer::trimReversePEAdapter(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    uint32_t cycles,
    FragmentMetadata & reverseFragment,
    Cigar &cigarBuffer) const
{
    if (cycles <= reverseFragment.leftClipped())
    {
        return false;
    }

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(reverseFragment.getCluster().getId(), " trim rv:" << reverseFragment << " cycles:" << cycles);
    const unsigned cigarOffset = cigarBuffer.size();
    uint64_t clip = 0;
    uint64_t align = 0;
    reference::ReferencePosition fwPos = reverseFragment.getFStrandReferencePosition();
    bool reverse = reverseFragment.reverse;

    Cigar::const_iterator current = reverseFragment.cigarBegin();
    for(;
        reverseFragment.cigarEnd() != current;
        ++current)
    {
        const Cigar::Component cigar = Cigar::decode(*current);
        if(Cigar::SOFT_CLIP == cigar.second)
        {
            clip += cigar.first;
            cycles -= std::min(cycles, cigar.first);
        }
        else if(Cigar::HARD_CLIP == cigar.second)
        {
            // Currently hard clips are only present in inversions. Too complicated to adjust, just bail out
            return false;
        }
        else if (Cigar::INSERT == cigar.second)
        {
            clip += cigar.first;
            cycles -= std::min(cycles, cigar.first);
        }
        else if (Cigar::DELETE == cigar.second)
        {
            fwPos += cigar.first;
        }
        else if (Cigar::BACK == cigar.second)
        {
            fwPos -= cigar.first;
        }
        else if (Cigar::FLIP == cigar.second || Cigar::CONTIG == cigar.second)
        {
            // too complicated to adjust, just bail out
            return false;
        }
        else if (Cigar::ALIGN == cigar.second)
        {
            if (cycles < cigar.first)
            {
                clip += cycles;
                align = cigar.first - cycles;
                fwPos += cycles;
                cycles = 0;
                break;
            }

            clip += cigar.first;
            fwPos += cigar.first;
            cycles -= cigar.first;
        }
        else
        {
            ISAAC_VERIFY_MSG(false, "Unexpected cigar operation:" << cigar);
        }
    }

    const Cigar::const_iterator oldCigarEnd = reverseFragment.cigarEnd();
    reverseFragment.resetAlignment();
    if (!align && oldCigarEnd == current)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(reverseFragment.getCluster().getId(), " fultrim:" << reverseFragment);
    }
    else
    {
        ISAAC_ASSERT_MSG(!cycles, "Some cycles left unclipped:" << cycles << " " << reverseFragment);

        if (clip)
        {
            cigarBuffer.addOperation(clip, Cigar::SOFT_CLIP);
        }
        if (align)
        {
            cigarBuffer.addOperation(align, Cigar::ALIGN);
        }
        cigarBuffer.addOperations(current + 1, oldCigarEnd);

        reverseFragment.updateAlignment(
            collectMismatchCycles_, alignmentCfg_, readMetadata, contigList,
            reverse, fwPos.getContigId(), fwPos.getPosition(),
            cigarBuffer, cigarOffset);
        reverseFragment.leftClipped() = std::max<unsigned short>(clip, reverseFragment.leftClipped());
        reverseFragment.incrementAdapterClip(clip);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(reverseFragment.getCluster().getId(), " trimmed:" << reverseFragment);
    }
    return true;
}

bool PeAdapterTrimmer::trimPEAdapterCycles(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const unsigned cycles,
    FragmentMetadata &fragment,
    Cigar &cigarBuffer) const
{
    return fragment.reverse ?
        trimReversePEAdapter(contigList, readMetadataList[fragment.getReadIndex()], cycles, fragment, cigarBuffer) :
        trimForwardPEAdapter(contigList, readMetadataList[fragment.getReadIndex()], cycles, fragment, cigarBuffer);
}


bool PeAdapterTrimmer::trimPEAdapter(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    FragmentMetadata & forwardFragment,
    FragmentMetadata & reverseFragment,
    Cigar &cigarBuffer) const
{
    bool ret = false;
    if (!forwardFragment.splitAlignment &&
        reverseFragment.getRStrandReferencePosition() >= forwardFragment.getFStrandReferencePosition() &&
        reverseFragment.getRStrandReferencePosition() < forwardFragment.getRStrandReferencePosition())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " trim fw:" << forwardFragment);
        forwardFragment = trimForwardPEAdapter(
            contigList, readMetadataList[forwardFragment.getReadIndex()], forwardFragment, reverseFragment.getRStrandReferencePosition(), cigarBuffer);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(forwardFragment.getCluster().getId(), " trimmed:" << forwardFragment);
        ret = true;
    }

    if (!reverseFragment.splitAlignment &&
        forwardFragment.getFStrandReferencePosition() > reverseFragment.getFStrandReferencePosition() &&
        forwardFragment.getFStrandReferencePosition() <= reverseFragment.getRStrandReferencePosition())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(reverseFragment.getCluster().getId(), " trim rv:" << reverseFragment);
        reverseFragment = trimReversePEAdapter(
            contigList, readMetadataList[reverseFragment.getReadIndex()], reverseFragment, forwardFragment.getFStrandReferencePosition(), cigarBuffer);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(reverseFragment.getCluster().getId(), " trimmed:" << reverseFragment);
        ret = true;
    }

    return ret;
}

bool PeAdapterTrimmer::checkTrimPEAdapter(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    FragmentMetadata &fragmentA,
    FragmentMetadata &fragmentB,
    Cigar &cigarBuffer) const
{
    if (trimPEAdapters_)
    {
        if (fragmentA.getContigId() == fragmentB.getContigId() &&
            !fragmentA.isSplit() && !fragmentB.isSplit())
        {
            if (fragmentA.isReverse() != fragmentB.isReverse())
            {
                if (!fragmentA.isReverse())
                {
                    return trimPEAdapter(contigList, readMetadataList, fragmentA, fragmentB, cigarBuffer);
                }
                else
                {
                    return trimPEAdapter(contigList, readMetadataList, fragmentB, fragmentA, cigarBuffer);
                }
            }
        }
    }
    return false;
}

void PeAdapterTrimmer::trimPEAdapterCycles(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const unsigned cycles,
    FragmentMetadataList &fragments,
    Cigar &cigarBuffer) const
{
    bool trimmed = false;
    for (FragmentMetadata &fragment : fragments)
    {
        trimmed |= trimPEAdapterCycles(contigList, readMetadataList, cycles, fragment, cigarBuffer);
    }

    if (trimmed)
    {
        fragments.erase(
            std::remove_if(fragments.begin(), fragments.end(), [](const FragmentMetadata &fragment){return !fragment.isAligned();}),
            fragments.end());
        std::sort(fragments.begin(), fragments.end());
        fragments.erase(std::unique(fragments.begin(), fragments.end()), fragments.end());
    }

}


} // namespace templateBuilder
} // namespace alignment
} // namespace isaac
