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
 ** \file ShadowAligner.cpp
 **
 ** \brief See ShadowAligned.hh
 ** 
 ** \author Come Raczy
 **/

#include <algorithm>
#include <boost/format.hpp>

#include "alignment/Quality.hh"
#include "alignment/templateBuilder/FragmentBuilder.hh"
#include "alignment/templateBuilder/ShadowAligner.hh"
#include "alignment/templateBuilder/FragmentSequencingAdapterClipper.hh"
#include "oligo/Kmer.hh"
#include "oligo/KmerGenerator.hpp"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

template <unsigned SHADOW_KMER_LENGTH>
ShadowAligner<SHADOW_KMER_LENGTH>::ShadowAligner(
    const bool collectMismatchCycles,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const bool noSmithWaterman,
    const bool splitAlignments,
    const AlignmentCfg &alignmentCfg,
    Cigar &cigarBuffer)
    : gappedMismatchesMax_(gappedMismatchesMax),
      smitWatermanGapsMax_(smitWatermanGapsMax),
      noSmithWaterman_(noSmithWaterman),
      splitAlignments_(splitAlignments),
      flowcellLayoutList_(flowcellLayoutList),
      ungappedAligner_(collectMismatchCycles, alignmentCfg),
      cigarBuffer_(cigarBuffer)
{
    static const std::size_t SHADOW_CANDIDATE_POSITIONS_MAX_EVER = 10000;
    shadowCandidatePositions_.reserve(SHADOW_CANDIDATE_POSITIONS_MAX_EVER);
}

template <unsigned SHADOW_KMER_LENGTH>
unsigned ShadowAligner<SHADOW_KMER_LENGTH>::hashShadowKmers(const std::vector<char> &sequence)
{
    shadowKmerPositions_.clear();
    // initialize all k-mers to the magic value -1 (NOT_FOUND)
    shadowKmerPositions_.resize(shadowKmerCount_, -1);
    // 
    oligo::KmerGenerator<SHADOW_KMER_LENGTH, unsigned, std::vector<char>::const_iterator> kmerGenerator(sequence.begin(), sequence.end());
    unsigned positionsCount = 0;
    unsigned kmer;
    std::vector<char>::const_iterator position;
    while (kmerGenerator.next(kmer, position))
    {
        if (-1 == shadowKmerPositions_[kmer])
        {
            shadowKmerPositions_[kmer] = (position - sequence.begin());
            ++positionsCount;
        }
    }
    return positionsCount;
}

template <unsigned SHADOW_KMER_LENGTH>
bool ShadowAligner<SHADOW_KMER_LENGTH>::findShadowCandidatePositions(
    const std::pair<int64_t, int64_t> &alignmentStartPositionRange,
    const int64_t referenceOffset,
    const reference::Contig::const_iterator referenceBegin,
    const reference::Contig::const_iterator referenceEnd,
    const std::vector<char> &shadowSequence,
    std::vector<int64_t> &shadowCandidatePositions)
{
    hashShadowKmers(shadowSequence);

    // find matching positions in the reference by k-mer comparison
    oligo::KmerGenerator<SHADOW_KMER_LENGTH, unsigned, reference::Contig::const_iterator> kmerGenerator(referenceBegin, referenceEnd);
    unsigned kmer;
    reference::Contig::const_iterator position = referenceBegin;
    while(kmerGenerator.next(kmer, position))
    {
        if (-1 != shadowKmerPositions_[kmer])
        {
            const int64_t candidatePosition = position - referenceBegin - shadowKmerPositions_[kmer] + referenceOffset;

            // avoid positions that will place mate outside the requested range. This can happen if rightmost k-mer of the mate matches the
            // first kmer of the reference like so:
            // <MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            //                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
            if (alignmentStartPositionRange.first <= candidatePosition && alignmentStartPositionRange.second >= candidatePosition)
            {
                // avoid spurious repetitions of start positions
                if (shadowCandidatePositions.empty() || shadowCandidatePositions.back() != candidatePosition)
                {
                    if (shadowCandidatePositions.size() == shadowCandidatePositions.capacity())
                    {
                        // too many candidate positions. Just stop here. The alignment score will be miserable anyway.
                        break;
                    }
                    shadowCandidatePositions.push_back(candidatePosition);
    //                ISAAC_THREAD_CERR << "kmer: " << oligo::Bases<2, oligo::Kmer>(kmer, shadowKmerLength_) << " " << candidatePosition << std::endl;
                }
            }
        }
    }

    // remove duplicate positions
    if (!shadowCandidatePositions.empty())
    {
        std::sort(shadowCandidatePositions.begin(), shadowCandidatePositions.end());
        shadowCandidatePositions.erase(std::unique(shadowCandidatePositions.begin(), shadowCandidatePositions.end()), shadowCandidatePositions.end());
        return true;
    }

    return false;
}

/**
 * \brief if the best template is longer than the dominant template, attempt to rescue shadow
 *        within a wider range
 * \param bestTemplateLength 0 indicates there is no best template
 */
std::pair <int64_t, int64_t> calculateShadowRescueRange(
    const FragmentMetadata &orphan,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    const Cluster &cluster = orphan.getCluster();

    const unsigned readLengths[] = {cluster[0].getLength(), cluster[1].getLength()};
    int64_t shadowMinPosition = templateLengthStatistics.mateMinPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths);
    int64_t shadowMaxPosition = templateLengthStatistics.mateMaxPosition(orphan.readIndex, orphan.reverse, orphan.position, readLengths);

    const std::pair <int64_t, int64_t> ret(shadowMinPosition, shadowMaxPosition);
    return ret;
}

template <unsigned SHADOW_KMER_LENGTH>
bool ShadowAligner<SHADOW_KMER_LENGTH>::findShadowCandidatePositions(
    const FragmentMetadata& orphan,
    const TemplateLengthStatistics& templateLengthStatistics,
    const Read& shadowRead, const bool shadowReverse,
    const reference::Contig& contig,
    std::vector<int64_t> &shadowCandidatePositions)
{
    const std::pair<int64_t, int64_t> shadowRescueRange = calculateShadowRescueRange(orphan, templateLengthStatistics);
    if (shadowRescueRange.second < shadowRescueRange.first)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(),
            "    Rescuing impossible: shadowMaxPosition < shadowMinPosition " << shadowRescueRange.second << " < " << shadowRescueRange.first);
        return false;
    }
    if (shadowRescueRange.second + 1 + shadowRead.getLength() < 0)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(),
            "    Rescuing impossible: shadowMaxPosition + 1 + shadowRead.getLength() < 0 " << (shadowRescueRange.second + 1 + shadowRead.getLength()) << "%l < 0");
        return false;
    }
    // find all the candidate positions for the shadow on the identified reference region
    shadowCandidatePositions.clear();
    const std::vector<char>& shadowSequence = shadowReverse ? shadowRead.getReverseSequence() : shadowRead.getForwardSequence();
    const int64_t candidatePositionOffset = std::min((int64_t) (contig.size()), std::max(int64_t(0), shadowRescueRange.first));
    findShadowCandidatePositions(
        shadowRescueRange,
        candidatePositionOffset,
        contig.begin() + candidatePositionOffset,
        contig.begin() + std::min((int64_t) (contig.size()), std::max(int64_t(0), shadowRescueRange.second) + 1),
        shadowSequence,
        shadowCandidatePositions);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
        orphan.getCluster().getId(),
        "findShadowCandidatePositions found " << shadowCandidatePositions.size() << " positions in range [" << (candidatePositionOffset) << ";" << (std::min((int64_t)contig.size(), shadowRescueRange.second + 1)) << "]");
    return !shadowCandidatePositions.empty();
}

template <unsigned SHADOW_KMER_LENGTH>
bool ShadowAligner<SHADOW_KMER_LENGTH>::alignCandidates(
    const FragmentMetadata& orphan,
    const bool reverse,
    const reference::ContigList& contigList,
    const flowcell::ReadMetadata& readMetadata,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    FragmentMetadataList& shadowList)
{
    // align the shadow to the candidate positions and keep the best fragment
    shadowList.clear();
    for (int64_t strandPosition : shadowCandidatePositions_)
    {
        FragmentMetadata fragment(
            &orphan.getCluster(), readMetadata.getIndex(), 0, 0,
            reverse, orphan.contigId, strandPosition, contigList[orphan.contigId].isDecoy());
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing: " << fragment);
        adapterClipper.checkInitStrand(fragment, contigList[orphan.contigId]);
        if (ungappedAligner_.alignUngapped(fragment, cigarBuffer_, readMetadata, adapterClipper, contigList))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing: Aligned: " << fragment);
            shadowList.push_back(fragment);
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing: Unaligned: " << fragment);
        }
    }

    return !shadowList.empty();
}

/**
 * \return false when no reasonable placement for shadow found. If at that point the shadowList is not empty,
 * this means that the shadow falls at a repetitive region and rescuing should not be considered
 */
template <unsigned SHADOW_KMER_LENGTH>
bool ShadowAligner<SHADOW_KMER_LENGTH>::rescueShadows(
    const reference::ContigList &contigList,
    const FragmentMetadata &orphan,
    const unsigned shadowsMax,
    FragmentMetadataList &shadowList,
    const flowcell::ReadMetadata &shadowReadMetadata,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const TemplateLengthStatistics &templateLengthStatistics)
{
    if (!templateLengthStatistics.isCoherent())
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescuing impossible. Incoherent tls");
        return false;
    }
    const Cluster &cluster = orphan.getCluster();
    const unsigned shadowReadIndex = (orphan.readIndex + 1) % 2;
    const Read &shadowRead = cluster[shadowReadIndex];
    // identify the orientation and range of reference positions of the orphan
    const bool reverse = templateLengthStatistics.mateOrientation(orphan.readIndex, orphan.reverse);
    if (!findShadowCandidatePositions(orphan, templateLengthStatistics, shadowRead, reverse, contigList[orphan.contigId], shadowCandidatePositions_) ||
        shadowsMax < shadowCandidatePositions_.size())
    {
        return false;
    }

    // align the shadow to the candidate positions and keep the best fragment
    if (!alignCandidates(orphan, reverse, contigList, shadowReadMetadata, adapterClipper, shadowList))
    {
        return false;
    }

    putBestOnTop<false>(shadowList);
    return true;
}

template class ShadowAligner<7>;
template class ShadowAligner<8>;

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac
