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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentBuilder.hh
 ** 
 ** \author Come Raczy
 **/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "common/Exceptions.hh"
#include "alignment/Mismatch.hh"
#include "alignment/templateBuilder/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

//std::array<std::atomic<std::size_t>, 7> FragmentBuilder::seedMatchCounts_;
//bool FragmentBuilder::countsTraced_ = false;

FragmentBuilder::FragmentBuilder(
    const bool collectMismatchCycles,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const unsigned repeatThreshold,
    const unsigned seedLength,
    const unsigned maxSeedsPerMatch,
    const unsigned seedRepeatThreshold,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const unsigned smithWatermanGapSizeMax,
    const bool noSmithWaterman,
    const bool splitAlignments,
    const AlignmentCfg &alignmentCfg,
    Cigar &cigarBuffer,
    const bool reserveBuffers)
    : repeatThreshold_(repeatThreshold)
    , gappedMismatchesMax_(gappedMismatchesMax)
    , smitWatermanGapsMax_(smitWatermanGapsMax)
    , noSmithWaterman_(noSmithWaterman)
    , smartSmithWaterman_(smartSmithWaterman)
    , splitAlignments_(splitAlignments)
    , alignmentCfg_(alignmentCfg)
    , cigarBuffer_(cigarBuffer)
    , ungappedAligner_(collectMismatchCycles, alignmentCfg_)
    , gappedAligner_(collectMismatchCycles, flowcellLayoutList, smartSmithWaterman, smithWatermanGapSizeMax, alignmentCfg_)
    , matchLists_(maxSeedsPerMatch + 1)
{
//    if (reserveBuffers)
    {
//        ISAAC_TRACE_STAT("FragmentBuilder before matchLists_.reserve with maxSeedsPerRead:" << maxSeedsPerRead);
        unsigned seedsPerMatch = 0;
        for (Matches &matches : matchLists_)
        {
            if (seedsPerMatch)
            {
                // max number of matches we get from a seed is seedRepeatThreshold. If the seed is combined with any other seed
                // the number of matches of a pair can only get smaller. So, use the number of non-overlapping tuples to high-bound the
                // matches coming out of seedsPerMatch
                matches.reserve(seedRepeatThreshold * maxSeedsPerMatch / seedsPerMatch);
            }
            ++seedsPerMatch;
        }
//        ISAAC_TRACE_STAT("FragmentBuilder before mergeBuffer_.reserve");
        // have buffer for one more than repeatThreshold_ as we have to insert new before we can remove the worst one
//        ISAAC_TRACE_STAT("FragmentBuilder before bestMatches_.reserve");
        bestMatches_.reserve(repeatThreshold_ + 1);
//        ISAAC_TRACE_STAT("FragmentBuilder after bestMatches_.reserve");
        ISAAC_ASSERT_MSG(MIN_CANDIDATES < repeatThreshold_, "repeatThreshold_ " << repeatThreshold_ << " is less than MIN_CANDIDATES " << MIN_CANDIDATES);
        oneFragmentCigarBuffer_.reserve(10240);

        fwMergeBuffers_.resize(matchLists_.size());
        rvMergeBuffers_.resize(matchLists_.size());

        const std::size_t maxSeedsPerRead = flowcell::getMaxReadLength(flowcellLayoutList) / seedLength;
        // [0] is being used for all sorts of temporary stuff. It should be able to hold max matches
        fwMergeBuffers_[0].reserve(seedRepeatThreshold * maxSeedsPerRead);
        rvMergeBuffers_[0].reserve(seedRepeatThreshold * maxSeedsPerRead);
        // [1] accumulates matches of all tried seeds. It should be able to hold max matches
        fwMergeBuffers_[1].reserve(seedRepeatThreshold * maxSeedsPerRead);
        rvMergeBuffers_[1].reserve(seedRepeatThreshold * maxSeedsPerRead);
        for (unsigned seedsPerMatch = 2; maxSeedsPerMatch >= seedsPerMatch; ++seedsPerMatch)
        {
            ISAAC_TRACE_STAT("before ClusterHashMatchFinder::ClusterHashMatchFinder seedsPerMatch:" << seedsPerMatch << "/" << maxSeedsPerMatch);
            // max number of matches we get from a seed is seedRepeatThreshold. If the seed is combined with any other seed
            // the number of matches of a pair can only get smaller. So, use the number of non-overlapping tuples to high-bound the
            // matches coming out of seedsPerMatch
            fwMergeBuffers_[seedsPerMatch].reserve(seedRepeatThreshold * maxSeedsPerMatch / seedsPerMatch);
            rvMergeBuffers_[seedsPerMatch].reserve(seedRepeatThreshold * maxSeedsPerMatch / seedsPerMatch);
            ISAAC_TRACE_STAT("ClusterHashMatchFinder::ClusterHashMatchFinder seedsPerMatch:" << seedsPerMatch << "/" << maxSeedsPerMatch);
        }
    }

}

//unsigned iSAAC_PROFILING_NOINLINE countMismatches(
//    const reference::Contig &contig,
//    const Read &read,
//    const bool reverse,
//    int64_t alignmentPosition)
//{
//    const std::vector<char> &sequence = read.getStrandSequence(reverse);
////    const std::vector<char> &quality = read.getStrandQuality(reverse);
//
//    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
//    std::vector<char>::const_iterator sequenceEnd = sequence.end();
//
//    templateBuilder::AlignerBase::clipReference(contig.size(), alignmentPosition, sequenceBegin, sequenceEnd);
////    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
//
//    return alignment::countMismatches(sequenceBegin, sequenceEnd, contig.begin() + alignmentPosition);
//}

void iSAAC_PROFILING_NOINLINE FragmentBuilder::updateBestMatches(
    const Match &match,
    const unsigned mismatches,
    std::vector<BestMatch> &bestMatches) const
{
    if (bestMatches.size() < (bestMatches.capacity() - 1) || mismatches < bestMatches.front().mismatches_)
    {
        bestMatches.push_back(BestMatch(match, mismatches));
        std::push_heap(bestMatches.begin(), bestMatches.end(),
                       [](const BestMatch &left, const BestMatch &right)
                       {return left.mismatches_ < right.mismatches_;});
        if (bestMatches.capacity() == bestMatches.size())
        {
            std::pop_heap(bestMatches.begin(), bestMatches.end(),
                           [](const BestMatch &left, const BestMatch &right)
                           {return left.mismatches_ < right.mismatches_;});
            bestMatches.pop_back();
        }
    }
}

unsigned countMismatches(
    const Match& match,
    const reference::ContigList& contigList,
    const Read& read)
{
    ISAAC_ASSERT_MSG(contigList.endOffset() >= match.contigListOffset_, "match.contigListOffset_ is outside valid range:" << match);
    const int64_t alignmentReferenceOffset = match.contigListOffset_;
    ISAAC_ASSERT_MSG(0 <= alignmentReferenceOffset, "alignmentPosition is negative:" << match);

    const std::vector<char> &sequence = read.getStrandSequence(match.reverse_);
    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    return alignment::countMismatches(sequenceBegin, sequenceEnd, contigList.referenceBegin() + alignmentReferenceOffset);
}
//
//void prefetch(
//    const Match& match,
//    const reference::ContigList& contigList)
//{
//    const int64_t alignmentReferenceOffset = match.contigListOffset_;
//    __builtin_prefetch(&*(contigList.referenceBegin() + alignmentReferenceOffset), 0, 0);
//    __builtin_prefetch(&*(contigList.referenceBegin() + alignmentReferenceOffset + 32), 0, 0);
//    __builtin_prefetch(&*(contigList.referenceBegin() + alignmentReferenceOffset + 64), 0, 0);
//    __builtin_prefetch(&*(contigList.referenceBegin() + alignmentReferenceOffset + 96), 0, 0);
//}

void FragmentBuilder::collectBestMatches(
    const Matches& matches,
    const reference::ContigList& contigList, const Read &read,
//    std::size_t &counts,
    std::vector<BestMatch>& bestMatches) const
{
//    ISAAC_ASSERT_MSG(matches.end() == std::adjacent_find(matches.begin(), matches.end()), "Duplicate matches unexpected:" << *std::adjacent_find(matches.begin(), matches.end()));

//    const int PREFETCH_OFFSET = 2;
    for (Matches::const_iterator it = matches.begin(); matches.end() != it; ++it)
    {
//        if (std::distance(it, matches.end()) > PREFETCH_OFFSET)
//        {
//            prefetch(*(it + PREFETCH_OFFSET), contigList);
//        }
        const Match& match = *it;
//        ISAAC_THREAD_CERR << match << std::endl;
        const unsigned mismatches = countMismatches(match, contigList, read);

//        ++counts;
        updateBestMatches(match, mismatches, bestMatches);
    }
}

//void FragmentBuilder::updateHitStats(const common::StaticVector<int, CHECK_MATCH_GROUPS_MAX>& seedCounts,
//                                     std::size_t counts[CHECK_MATCH_GROUPS_MAX]) const
//{
//    if (seedCounts.size() == 2)
//    {
//        matchSeedCounts_[seedCounts[0]][seedCounts[1]].first += counts[0];
//        matchSeedCounts_[seedCounts[0]][seedCounts[1]].second += counts[1];
//    }
//    else
//    {
//        ISAAC_ASSERT_MSG(1 == seedCounts.size(), "Empty seedCounts");
//        matchSeedCounts_[seedCounts[0]][0].first += counts[0];
//    }
//}

/**
 * \return true if at least one fragment was built.
 */
AlignmentType FragmentBuilder::findBestMatches(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const Cluster &cluster,
    const MatchLists &matchLists,
    std::vector<BestMatch> &bestMatches) const
{
    bestMatches.clear();
//    common::StaticVector<int, CHECK_MATCH_GROUPS_MAX> seedCounts;
//    std::size_t counts[CHECK_MATCH_GROUPS_MAX] = {0, 0};
//    std::size_t countsIndex = 0;

    for (unsigned supportingSeeds = matchLists.size() - 1; supportingSeeds; --supportingSeeds)
    {
        const Matches &matches = matchLists[supportingSeeds];
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    supportingSeeds:" << supportingSeeds << " matches:" << matches.size());
//        if (!matches.empty() && !countsIndex)
//        {
//            countsIndex = supportingSeeds;
//        }
//        seedMatchCounts_[countsIndex] += matches.size();
//        if (!matches.empty())
//        {
//            seedCounts.push_back(supportingSeeds);
//        }
//        // this may reduce number of random memory page requests if several positions hit the same page
//        std::sort(matches.begin(), matches.end(), [](const Match &left, const Match &right)
//                  {return left.contigListOffset_ < right.contigListOffset_;});

        collectBestMatches(
            matches, contigList, cluster.at(readMetadata.getIndex()),
            //counts[countsIndex],
            bestMatches);
//        countsIndex += !matches.empty();
    }

//    updateHitStats(seedCounts, counts);

    return bestMatches.empty() ? Nm : Normal;
}

/**
 * \return true if at least one fragment was built.
 */
FragmentMetadata  FragmentBuilder::makeAlignment(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const Cluster &cluster,
    FragmentSequencingAdapterClipper &adapterClipper,
    const Match &match,
    const unsigned repeatSeeds,
    Cigar &cigarBuffer) const
{
    const unsigned contigId = contigList.contigIdFromOffset(match.contigListOffset_);
    FragmentMetadata fragmentMetadata(
        &cluster, readMetadata.getIndex(), 0, 0, match.reverse_,
        contigId,
        contigList.positionFromOffset(match.contigListOffset_),
        contigList[contigId].isDecoy(),
        repeatSeeds);

    adapterClipper.checkInitStrand(fragmentMetadata, contigList.at(fragmentMetadata.contigId));
    ungappedAligner_.alignUngapped(fragmentMetadata, cigarBuffer, readMetadata, adapterClipper, contigList);
    return fragmentMetadata;
}

/**
 * \return true if at least one fragment was built.
 */
bool  iSAAC_PROFILING_NOINLINE FragmentBuilder::makeBestAlignments(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const Cluster &cluster,
    const bool withGaps,
    const std::vector<BestMatch> &bestMatches,
    const unsigned repeatSeeds,
    FragmentSequencingAdapterClipper &adapterClipper,
    FragmentMetadataList &fragments) const
{
    std::size_t secondBestMismatches = 0;
    bool perfectFound = false;
    for (const BestMatch &bestMatch : bestMatches)
    {
        if (!secondBestMismatches)
        {
            if (bestMatch.mismatches_ > bestMatches.front().mismatches_)
            {
                secondBestMismatches = bestMatch.mismatches_;
            }
        }
        else if (bestMatch.mismatches_ > secondBestMismatches)
        {
            // keep best and second best only. The difference between the two (apart from repeats) is what ultimately defines alignment score
            break;
        }

        const FragmentMetadata fragment = makeAlignment(contigList, readMetadata, cluster, adapterClipper, bestMatch.match_, repeatSeeds, cigarBuffer_);
        if (fragment.isAligned())
        {
            fragments.push_back(fragment);
            perfectFound |= !fragments.back().mismatchCount;
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    FragmentBuilder::makeBestAlignments " << contigList[fragments.back().contigId] << " " << bestMatch.match_<< " " << fragments.back());
        }

//        ++alignedSeedCounts_[bestMatch.match_.seedCount_];
    }

    bool gappedFound = false;
    // having perfect alignments means no need to spend time on trying to improve the imperfect ones.
    if (!noSmithWaterman_ && withGaps && (!perfectFound || !smartSmithWaterman_))
    {
        // If there are still bad alignments, try to do expensive smith-waterman on them.
        if (gappedAligner_.realignBadUngappedAlignments(
            gappedMismatchesMax_, smitWatermanGapsMax_, contigList, readMetadata, fragments, adapterClipper, cigarBuffer_))
        {
            const FragmentMetadataList::iterator bestGapped = std::min_element(fragments.begin(), fragments.end(), FragmentMetadata::bestGappedLess);
            const FragmentMetadataList::const_iterator secondBestGapped = std::min_element(bestGapped + 1, fragments.end(), FragmentMetadata::bestGappedLess);
            if (!bestGapped->gapCount && (fragments.end() == secondBestGapped || !secondBestGapped->gapCount))
            {
                fragments.erase(std::remove_if(fragments.begin(), fragments.end(), [](const FragmentMetadata &fragment){return fragment.gapCount;}), fragments.end());
            }
            else
            {
                gappedFound = true;
            }
        }
    }

    std::sort(fragments.begin(), fragments.end(), gappedFound ? FragmentMetadata::bestGappedLess : FragmentMetadata::bestUngappedLess);

    unsigned bestCount = 0;
    for (const FragmentMetadata &fragment : fragments)
    {
        if (FragmentMetadata::alignmentsEquivalent(fragment, fragments.front()))
        {
            ++bestCount;
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    FragmentBuilder::makeBestAlignments bestCount:" << bestCount << " " << fragment);
        }
        else
        {
            break;
        }
    }

    for (FragmentMetadata &fragment : fragments)
    {
        if (FragmentMetadata::alignmentsEquivalent(fragment, fragments.front()))
        {
            fragment.repeatCount = bestCount;
        }
        else
        {
            break;
        }
    }

    return !fragments.empty();
}

/**
 * \return true if at least one fragment was built.
 */
AlignmentType FragmentBuilder::findBestAlignments(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const Cluster &cluster,
    const bool withGaps,
    const MatchLists &matchLists,
    const unsigned repeatSeeds,
    FragmentMetadataList &fragments) const
{
    const AlignmentType res = findBestMatches(contigList, readMetadata, cluster, matchLists, bestMatches_);
    if (Normal != res)
    {
        return res;
    }

//    std::sort(bestMatches_.begin(), bestMatches_.end(), [](const BestMatch &left, const BestMatch &right)
//              {return left.mismatches_ < right.mismatches_ || (left.mismatches_ == right.mismatches_ && left.match_.seedCount_ > right.match_.seedCount_);});
    if (makeBestAlignments(contigList, readMetadata, cluster, withGaps, bestMatches_, repeatSeeds, adapterClipper, fragments))
    {
//        ISAAC_ASSERT_MSG("FC:28451918" != std::string(cluster.nameBegin(), cluster.nameEnd()), cluster);

        return Normal;
    }

    return Nm;
}


//std::array<std::array<std::pair<std::atomic<std::size_t>, std::atomic<std::size_t> >, 7>, 7> FragmentBuilder::matchSeedCounts_;
//std::array<std::atomic<std::size_t>, 7> FragmentBuilder::alignedSeedCounts_;
//bool FragmentBuilder::countsTraced_ = false;

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac
