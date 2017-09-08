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
 ** \file UngappedAligner.cpp
 **
 ** \brief See UngappedAligner.hh
 ** 
 ** \author Roman Petrovski
 **/
#include "alignment/templateBuilder/UngappedAligner.hh"
#include "alignment/Mismatch.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

UngappedAligner::UngappedAligner(
    const bool collectMismatchCycles,
    const AlignmentCfg &alignmentCfg)
    : AlignerBase(collectMismatchCycles, alignmentCfg)
{
}

unsigned UngappedAligner::alignUngapped(
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer,
    const flowcell::ReadMetadata &readMetadata,
    const templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::ContigList &contigList
    ) const
{
    const unsigned cigarOffset = cigarBuffer.size();

// Don't reset alignment to preserve the seed-based anchors.
//    fragmentMetadata.resetAlignment();
    ISAAC_ASSERT_MSG(!fragmentMetadata.isAligned(), "alignUngapped is expected to be performed on a clean fragment");
    fragmentMetadata.resetClipping();
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignUngapped: after resetClipping: " << fragmentMetadata);

    const reference::Contig &contig = contigList[fragmentMetadata.contigId];

    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const reference::Contig &reference = contig;

    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    adapterClipper.clip(contig, fragmentMetadata, sequenceBegin, sequenceEnd);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignUngapped: after adapterClipper.clip: " << fragmentMetadata);
    clipReadMasking(read, fragmentMetadata, sequenceBegin, sequenceEnd);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignUngapped: after clipReadMasking: " << fragmentMetadata);

    clipReference(reference.size(), fragmentMetadata.position, sequenceBegin, sequenceEnd);

    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
    if (firstMappedBaseOffset)
    {
        cigarBuffer.addOperation(firstMappedBaseOffset, Cigar::SOFT_CLIP);
    }

    const unsigned mappedBases = std::distance(sequenceBegin, sequenceEnd);
    if (mappedBases)
    {
        const Cigar::OpCode opCode = Cigar::ALIGN;
        cigarBuffer.addOperation(mappedBases, opCode);
    }

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
    }

    const unsigned ret = updateFragmentCigar(
        readMetadata, contigList, //contigAnnotations,
        fragmentMetadata,
        fragmentMetadata.reverse, fragmentMetadata.contigId, fragmentMetadata.position, cigarBuffer, cigarOffset);

    if (!ret)
    {
        fragmentMetadata.setUnaligned();
    }

    return ret;
}

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac
