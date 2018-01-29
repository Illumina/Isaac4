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
 ** \file TemplateBuilder.cpp
 **
 ** \brief See TemplateBuilder.hh
 **
 ** \author Roman Petrovski
 **/

#include <numeric>
#include <limits>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/Mismatch.hh"
#include "alignment/templateBuilder/FragmentBuilder.hh"
#include "alignment/Quality.hh"
#include "alignment/TemplateBuilder.hh"
#include "oligo/KmerGenerator.hpp"

namespace isaac
{
namespace alignment
{

using namespace templateBuilder;

//PairClassCounts TemplateBuilder::pairClassCounts_;

TemplateBuilder::TemplateBuilder(
    const bool collectMismatchCycles,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const std::size_t candidateMatchesMax,
    const unsigned repeatThreshold,
    const unsigned seedLength,
    const unsigned maxSeedsPerMatch,
    const unsigned matchFinderTooManyRepeats,
    const unsigned matchFinderWayTooManyRepeats,
    const unsigned matchFinderShadowSplitRepeats,
    const bool scatterRepeats,
    const bool rescueShadows,
    const bool trimPEAdapters,
    const bool anchorMate,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const bool smartSmithWaterman,
    const unsigned smithWatermanGapSizeMax,
    const bool splitAlignments,
    const AlignmentCfg &alignmentCfg,
    const DodgyAlignmentScore dodgyAlignmentScore,
    const unsigned anomalousPairHandicap,
    const bool reserveBuffers)
    : repeatThreshold_(repeatThreshold)
    , seedLength_(seedLength)
    , matchFinderTooManyRepeats_(matchFinderTooManyRepeats)
    , matchFinderWayTooManyRepeats_(matchFinderWayTooManyRepeats)
    , matchFinderShadowSplitRepeats_(matchFinderShadowSplitRepeats)
    , scatterRepeats_(scatterRepeats)
    , rescueShadows_(rescueShadows)
    , anchorMate_(anchorMate)
    , dodgyAlignmentScore_(dodgyAlignmentScore)
    , anomalousPairHandicap_(pow10(double(anomalousPairHandicap) / 10.0) - 1.0)//anomalousPairHandicap)
    , flowcellLayoutList_(flowcellLayoutList)
    , smitWatermanGapsMax_(smitWatermanGapsMax)
    , splitAlignments_(splitAlignments)
    , alignmentCfg_(alignmentCfg)
    , peAdapterTrimmer_(collectMismatchCycles, trimPEAdapters, alignmentCfg_)
    , cigarBuffer_()
    , fragmentBuilder_(
        collectMismatchCycles, flowcellLayoutList, repeatThreshold_, seedLength_, maxSeedsPerMatch,
        std::max(matchFinderTooManyRepeats, std::max(matchFinderWayTooManyRepeats, matchFinderShadowSplitRepeats)),
        gappedMismatchesMax, smitWatermanGapsMax,
        smartSmithWaterman, smithWatermanGapSizeMax, !smitWatermanGapsMax_, splitAlignments,
        alignmentCfg_, cigarBuffer_, reserveBuffers)
    , shadowAligner_(collectMismatchCycles, flowcellLayoutList,
                     gappedMismatchesMax, smitWatermanGapsMax, smartSmithWaterman, !smitWatermanGapsMax_, splitAlignments, alignmentCfg_, cigarBuffer_)
    , splitReadAligner_(collectMismatchCycles, alignmentCfg_)
    , bestCombinationPairInfo_(0)
    , bestRescuedPair_(0)

{
//    if (reserveBuffers)
    {
//        ISAAC_TRACE_STAT("TemplateBuilder before cigarBuffer_.reserve");
        // if smith waterman is enabled, each candidate may produce up to one more gapped cigar
        const bool gapsPossible = std::max<unsigned>(smitWatermanGapsMax, splitAlignments);
        cigarBuffer_.reserve(
            // Cigar::getMaxOperationsForReads already accounts for read pairing
            Cigar::getMaxOperationsForReads(flowcellLayoutList, false) +
            gapsPossible * Cigar::getMaxOperationsForReads(flowcellLayoutList, true) * repeatThreshold_ +
            rescueShadows_ ?
                Cigar::getMaxOperationsForReads(flowcellLayoutList, gapsPossible) *
                std::max(matchFinderTooManyRepeats, std::max(matchFinderWayTooManyRepeats, matchFinderShadowSplitRepeats)) : 0);

        ISAAC_TRACE_STAT("TemplateBuilder before shadowList_.reserve");
        // each shadow can have up to 1 gapped alignment + have some room left to mix seed and rescued candidates for scoring
        shadowList_[0].reserve(BEST_SHADOWS_TO_KEEP * 2 + TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING);
        shadowList_[1].reserve(BEST_SHADOWS_TO_KEEP * 2 + TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING);
        ISAAC_TRACE_STAT("TemplateBuilder before bestCombinationPairInfo_.reserve");
        bestCombinationPairInfo_.reserve(repeatThreshold_, 0);
        ISAAC_TRACE_STAT("TemplateBuilder before bestRescuedPair_.reserve");
        bestRescuedPair_.reserve(repeatThreshold_, repeatThreshold_);
        ISAAC_TRACE_STAT("TemplateBuilder before candidates_.reserve");
        // each candidate can have up to 1 gapped alignment
        std::for_each(
            candidates_.begin(), candidates_.end(), boost::bind(&FragmentMetadataList::reserve, _1,
            repeatThreshold_ + gapsPossible * repeatThreshold_));
        ISAAC_TRACE_STAT("TemplateBuilder after candidates_.reserve");
    }
}

/**
 * \return false if template needs to go into unaligned bin
 */
templateBuilder::AlignmentType TemplateBuilder::flagDodgyTemplate(BamTemplate &bamTemplate) const
{
    ISAAC_ASSERT_MSG(2 == bamTemplate.getFragmentCount(), "Only paired can be flagged as dodgy");
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getFragmentMetadata(0).getCluster().getId(), "flagDodgyTemplate with " << dodgyAlignmentScore_ << " : " << bamTemplate);

    bamTemplate.resetAlignmentScore();
    if (DODGY_ALIGNMENT_SCORE_UNALIGNED == dodgyAlignmentScore_)
    {
        // both must sort into unaligned bin. setUnaligned will not do it.
        bamTemplate.setNoMatch();
        return Rm;
    }
    else
    {
        bamTemplate.setDodgy(dodgyAlignmentScore_);
    }
    return Normal;
}

templateBuilder::AlignmentType TemplateBuilder::buildCombinationTemplate(
    const reference::ContigList &contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList &readMetadataList,
    FragmentMetadataLists &fragments,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics,
    BamTemplate &bamTemplate) const
{
    // Initialize the bamTemplate_ with unaligned fragments for the cluster
    bamTemplate.reset(readMetadataList, cluster);
    if (2 == readMetadataList.size())
    {
        ISAAC_ASSERT_MSG(!fragments[0].empty() || !fragments[1].empty(), "Both paired fragment lists are empty");
        ISAAC_ASSERT_MSG(2 == readMetadataList.size(), "Paired fragment lists without paired readMetadata");
        if (!fragments[0].empty() && !fragments[1].empty())
        {
            if (!pickBestPair(contigList, readMetadataList, rog, fragments, templateLengthStatistics, bamTemplate))
            {
                return templateBuilder::Rm;
            }
        }
        else if (fragments[0].empty() != fragments[1].empty())
        {
            const unsigned orphanIndex = fragments[0].empty();
            buildSingletonShadowTemplate(
                rog, templateLengthStatistics, fragments[orphanIndex], *getBestFragment(fragments[orphanIndex]), orphanIndex, bamTemplate);
        }
    }
    else if(1 == readMetadataList.size() || 1 == fragments.size())
    {
        ISAAC_ASSERT_MSG(fragments[1].empty(), "With single-ended data expecting the fragment to be placed at index 0");
        if (!fragments[0].empty())
        {
            pickBestFragment(rog, fragments[0], bamTemplate);
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Single-ended fragment list is empty");
        }
    }
    else
    {
        ISAAC_ASSERT_MSG(false, (boost::format("TemplateBuilder supports at most 2 reads: %d reads found") % fragments.size()).str().c_str());
    }
    return templateBuilder::Normal;
}


FragmentMetadataList::const_iterator TemplateBuilder::getBestFragment(const FragmentMetadataList &fragments) const
{
    ISAAC_ASSERT_MSG(!fragments.empty(), "Empty list not allowed");

    const FragmentIterator bestAlignment = fragments.begin();
    const FragmentIterator notBestAlignment =
        std::find_if(bestAlignment + 1, fragments.end(),
                     [bestAlignment](const FragmentMetadata &fragment){return !FragmentMetadata::alignmentsEquivalent(*bestAlignment, fragment);});

    ISAAC_ASSERT_MSG(
        fragments.end() == notBestAlignment ||
        !notBestAlignment->isBetterGapped(*bestAlignment) || !notBestAlignment->isBetterUngapped(*bestAlignment),
        "Incorrect candidate order. Expected sorted best on top " << *bestAlignment<< " vs " << *notBestAlignment);

    const unsigned clusterId = bestAlignment->getCluster().getId();
    const unsigned repeatIndex = scatterRepeats_ ? (clusterId % std::distance(bestAlignment, notBestAlignment)) : 0;

    ISAAC_THREAD_CERR_DEV_TRACE("TemplateBuilder::getBestFragment returning repeat " <<
                                repeatIndex << " out of " << std::distance(bestAlignment, notBestAlignment) << " cluster id: " << clusterId << " " <<
                                *(bestAlignment + repeatIndex));

    return bestAlignment + repeatIndex;

}


/**
 * \return false if too many best pairs
 */
bool TemplateBuilder::locateBestAnchoredPair(
    const reference::ContigList &contigList,
    const RestOfGenomeCorrection &rog,
    const FragmentMetadataLists &fragments,
    const TemplateLengthStatistics &tls,
    unsigned bestPairsMax,
    BestPairInfo &ret) const
{
    ret.clear();
    unsigned bestPairsLeft = bestPairsMax;

//    const auto bestR1 = fragments[0].begin();
//    const auto bestR2 = fragments[1].begin();
    for(FragmentIterator r1Fragment = fragments[0].begin(); bestPairsLeft && fragments[0].end() != r1Fragment; ++r1Fragment)
    {
        //ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(r1Fragment->getCluster().getId(), " locateBestAnchoredPair r1Fragment:" << *r1Fragment);
//        ++pairClassCounts_.totalR1Candidates_;
//        if (r1Fragment->decoyAlignment)
//        {
//            ++pairClassCounts_.decoyR1Candidates_;
//        }
//        if (!FragmentMetadata::alignmentsEquivalent(*bestR1, *r1Fragment))
//        {
//            ISAAC_ASSERT_MSG(!r1Fragment->isBetterGapped(*bestR1) || !r1Fragment->isBetterUngapped(*bestR1), "Incorrect candidate order. Expected sorted best on top " << *bestR1 << " vs " << *r1Fragment);
//            break;
//        }
        for(FragmentIterator r2Fragment = fragments[1].begin(); bestPairsLeft && fragments[1].end() != r2Fragment; ++r2Fragment)
        {
            //ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(r2Fragment->getCluster().getId(), " locateBestAnchoredPair r2Fragment:" << *r2Fragment);
//            if (!FragmentMetadata::alignmentsEquivalent(*bestR2, *r2Fragment))
//            {
//                ISAAC_ASSERT_MSG(!r2Fragment->isBetterGapped(*bestR2) || !r2Fragment->isBetterUngapped(*bestR2), "Incorrect candidate order. Expected sorted best on top " << *bestR2<< " vs " << *r2Fragment);
//                break;
//            }
            int res = updateBestAnchoredPair(anomalousPairHandicap_, rog, tls, *r1Fragment, *r2Fragment, ret);
            // not strictly counting best, as there is no guarantee that they come in best to worst order, but should be
            // a reasonable approximation to avoid accumulating insane number of repeats.
            if (-1 == res)
            {
                bestPairsLeft = bestPairsMax - 1;
            }
            else if (res)
            {
                --bestPairsLeft;
            }
        }
    }

    if (!bestPairsLeft)
    {
        return false;
    }

    for (const FragmentMetadata &r1Fragment :fragments[0])
    {
        ret.appendSingleProbability(r1Fragment);
    }
    for (const FragmentMetadata &r2Fragment :fragments[1])
    {
        ret.appendSingleProbability(r2Fragment);
    }

    ISAAC_ASSERT_MSG(!ret.empty(), "Must pick one");
    return true;
}

/**
 * \brief   recomputes tail anchors and moves all semialigned to the end of shadowList.
 * @return  iterator to first semialigned or shadowList.end()
 */
FragmentMetadataList::iterator TemplateBuilder::pushSemialignedDown(
    const reference::ContigList& contigList,
    FragmentMetadataList& shadowList,
    std::size_t &semialignedHeadLength) const
{
    semialignedHeadLength = std::size_t(0) - 1;
    FragmentMetadataList::iterator ret = shadowList.end();
    for (FragmentMetadataList::iterator it = shadowList.begin(); ret != it;)
    {
        FragmentMetadata &shadow = *it;
// interesting idea but it anchors heads of reads that for some reason begin with NNNN
//        if (shadow.recomputeHeadAnchor<-1U, SEMIALIGNED_MATCHES_MIN>(contigList))
//        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadow.getCluster().getId(), "    UHA: " << shadow);
//            // if it is anchored at the head, it is very unlikely to be semialigned
//            ++it;
//            continue;
//        }
        // don't mess with gapped or clipped as they are difficult and not useful
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadow.getCluster().getId(), "    trying: " << shadow);
        if (!shadow.gapCount && !clippedByReference(contigList, shadow))
        {
            shadow.recomputeTailAnchor(contigList);

            if (SEMIALIGNED_MATCHES_MIN <= shadow.tailAnchor().length())
            {
                // find the end of assumedly clean tail so that we know the amount of cycles to use
                // for split alignment search
                const Anchor ha = shadow.computeAnchor<true, false, 0>(contigList, seedLength_);
                if (!ha.empty())
                {
                    semialignedHeadLength = std::min<std::size_t>(semialignedHeadLength,
                        shadow.isReverse() ?
                            shadow.getReadLength() - ha.second :
                            ha.first);
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadow.getCluster().getId(), "    UTA: " << shadow);
                    --ret;
                    std::swap(shadow, *ret);
                    continue;
                }
            }
        }
        ++it;
    }

    return ret;
}

void TemplateBuilder::buildSingletonShadowTemplate(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const TemplateLengthStatistics &templateLengthStatistics,
    const FragmentMetadataList &orphans,
    const FragmentMetadata &bestOrphan,
    const unsigned orphanIndex,
    BamTemplate &bamTemplate) const
{
    FragmentMetadata shadow = bestOrphan;
    // mark shadow as 'shadow' (singleton's position, etc)
    shadow.setUnaligned();
    shadow.readIndex = !orphanIndex;

    FragmentMetadata orphan = bestOrphan;

    orphan.alignmentScore = computeAlignmentScore(orphan, restOfGenomeCorrection, orphans);
    orphan.mapQ = alignmentScoreToMapq(orphan.alignmentScore);

    bamTemplate = BamTemplate(orphanIndex ? shadow : orphan, orphanIndex ? orphan : shadow, false);
}


void TemplateBuilder::pickRandomRepeatAlignment(
    const unsigned clusterId,
    const BestPairInfo &bestPair,
    BamTemplate &bamTemplate) const
{
    const unsigned repeatIndex = scatterRepeats_ ? clusterId % bestPair.repeatsCount() : 0;

//    const unsigned debugClass = bamTemplate.debugClass_;
    bamTemplate = bestPair.repeat(repeatIndex);
//    bamTemplate.debugClass_ = debugClass;

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(clusterId, "TemplateBuilder::pickRandomRepeatAlignment: Picked repeat " << repeatIndex <<
                                " out of " << bestPair.repeatsCount() << " " << bamTemplate);
}

void TemplateBuilder::scoreBestPair(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const FragmentMetadataLists &fragments,
    BestPairInfo &bestPair,
    BamTemplate &bamTemplate) const
{
    FragmentMetadata r1Alignment = bamTemplate.getFragmentMetadata(0);
    FragmentMetadata r2Alignment = bamTemplate.getFragmentMetadata(1);

    const double otherTemplateProbability =
        bestPair.sumUniquePairProbabilities(
            r1Alignment.logProbability + r2Alignment.logProbability,
            bestPair.repeatsCount(), bamTemplate.isProperPair());

    const unsigned pairScore = computeAlignmentScore(
        restOfGenomeCorrection.getRogCorrection(), bestPair.probability(), otherTemplateProbability);

    r1Alignment.alignmentScore = computeAlignmentScore(r1Alignment, restOfGenomeCorrection, fragments[0]);
    r2Alignment.alignmentScore = computeAlignmentScore(r2Alignment, restOfGenomeCorrection, fragments[1]);

    r1Alignment.mapQ = pickMapQ(
        r1Alignment.alignmentScore, r2Alignment.alignmentScore, bamTemplate.isProperPair(), pairScore);
    r2Alignment.mapQ = pickMapQ(
        r2Alignment.alignmentScore, r1Alignment.alignmentScore, bamTemplate.isProperPair(), pairScore);

    bamTemplate = BamTemplate(r1Alignment, r2Alignment, bamTemplate.isProperPair(), pairScore);

    ISAAC_ASSERT_MSG(
        bamTemplate.getAlignmentScore() < 4 || 1 == bestPair.repeatsCount() ,
        "alignment score too high for a repeat of " << bestPair.repeatsCount() << ":" << bamTemplate <<
        common::makeFastIoString(r1Alignment.getRead().getForwardSequence().begin(), r1Alignment.getRead().getForwardSequence().end()) <<
        "-" <<
        common::makeFastIoString(r2Alignment.getRead().getForwardSequence().begin(), r2Alignment.getRead().getForwardSequence().end()));
}

void TemplateBuilder::pickBestFragment(
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const FragmentMetadataList &fragmentList,
    BamTemplate &result) const
{
    ISAAC_ASSERT_MSG(!fragmentList.empty(), "Candidate list is empty");

    typedef FragmentMetadataList::const_iterator FragmentIterator;
    const FragmentIterator bestFragment = getBestFragment(fragmentList);
    FragmentMetadata fragment = *bestFragment;
    fragment.alignmentScore = computeAlignmentScore(fragment, restOfGenomeCorrection, fragmentList);
    fragment.mapQ = alignmentScoreToMapq(fragment.alignmentScore);
    result = BamTemplate(fragment);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentList.front().getCluster().getId(), "Single-ended template: " << result.getFragmentMetadata(0));
}

bool TemplateBuilder::pickBestPair(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    const RestOfGenomeCorrection &rog,
    FragmentMetadataLists &fragments,
    const TemplateLengthStatistics &templateLengthStatistics,
    BamTemplate &bamTemplate) const
{
    ISAAC_ASSERT_MSG(READS_IN_A_PAIR == fragments.size(), "TemplateBuilder::pickBestPair must be called for paired templates only");
    ISAAC_ASSERT_MSG(!fragments[0].empty() && !fragments[1].empty(), "TemplateBuilder::pickBestPair must be called for paired templates only");

    if (!locateBestAnchoredPair(contigList, rog, fragments, templateLengthStatistics, repeatThreshold_, bestCombinationPairInfo_))
    {
        return false;
    }

    bestCombinationPairInfo_.removeRepeatDuplicates();
    pickRandomRepeatAlignment(fragments[0][0].getCluster().getId(), bestCombinationPairInfo_, bamTemplate);

    FragmentMetadata read0 = bamTemplate.getFragmentMetadata(0);
    FragmentMetadata read1 = bamTemplate.getFragmentMetadata(1);
    if (peAdapterTrimmer_.checkTrimPEAdapter(
        contigList, readMetadataList, read0, read1, cigarBuffer_))
    {
        // If best choice ended up to be shorter than read length PE template, we need to trim adapter ends and recompute log probabilities
        peAdapterTrimmer_.trimPEAdapterCycles(
            contigList, readMetadataList, read0.highClipped, fragments[0], cigarBuffer_);
        peAdapterTrimmer_.trimPEAdapterCycles(
            contigList, readMetadataList, read1.highClipped, fragments[1], cigarBuffer_);

        ISAAC_ASSERT_MSG(!fragments[0].empty(), "All fragments are gone after cycle trimming");
        ISAAC_ASSERT_MSG(!fragments[1].empty(), "All fragments are gone after cycle trimming");

        std::sort(fragments[0].begin(), fragments[0].end(), FragmentMetadata::bestUngappedLess);
        std::sort(fragments[1].begin(), fragments[1].end(), FragmentMetadata::bestUngappedLess);

        if (!locateBestAnchoredPair(contigList, rog, fragments, templateLengthStatistics, repeatThreshold_, bestCombinationPairInfo_))
        {
            return false;
        }

        bestCombinationPairInfo_.removeRepeatDuplicates();
        pickRandomRepeatAlignment(fragments[0][0].getCluster().getId(), bestCombinationPairInfo_, bamTemplate);
    }

    scoreBestPair(rog, fragments, bestCombinationPairInfo_, bamTemplate);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragments.front().front().getCluster().getId(), "Pair-end  template well anchored: " << bamTemplate);
    return true;
}

void TemplateBuilder::trimShadowPairPEAdapaters(
    const reference::ContigList& contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList& readMetadataList,
    const TemplateLengthStatistics& tls,
    FragmentMetadata &orphan,
    FragmentMetadata &shadow,
    FragmentMetadataList &orphanList,
    FragmentMetadataList &shadowList,
    templateBuilder::BestPairInfo &bestRescuedPair) const
{
    if (peAdapterTrimmer_.checkTrimPEAdapter(contigList, readMetadataList, orphan, shadow, cigarBuffer_))
    {
        // If best choice ended up to be shorter than read length PE template, we need to trim adapter ends and recompute log probabilities
        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, orphan.highClipped, orphanList, cigarBuffer_);
        // recompute orphan scores after adapter trimming
        orphan.alignmentScore = computeAlignmentScore(orphan, rog, orphanList);
        orphan.mapQ = alignmentScoreToMapq(orphan.alignmentScore);

        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, shadow.highClipped, shadowList, cigarBuffer_);
        //repeat all pair building on trimmed fragments, but don't base best pair selection on the results
        //TODO: this is somewhat inefficient as we don't need to accumulate resulting pairs, only the pair probabilities
        ISAAC_ASSERT_MSG(!shadowList.empty(), "All fragments are gone after cycle trimming");
        buildRescuedPairs(anomalousPairHandicap_, contigList, rog, orphan, shadowList, tls, bestRescuedPair);
    }
}

void TemplateBuilder::trimDoublecheckedPairPEAdapaters(
    const reference::ContigList& contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList& readMetadataList,
    const TemplateLengthStatistics& tls,
    FragmentMetadata &orphan,
    FragmentMetadata &shadow,
    FragmentMetadata r0Orphan,
    FragmentMetadata r1Orphan,
    FragmentMetadataLists &orphanList,
    FragmentMetadataLists &shadowList,
    templateBuilder::BestPairInfo &bestRescuedPair) const
{
    if (peAdapterTrimmer_.checkTrimPEAdapter(contigList, readMetadataList, orphan, shadow, cigarBuffer_))
    {
        // If best choice ended up to be shorter than read length PE template, we need to trim adapter ends and recompute log probabilities
        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, orphan.highClipped, orphanList[orphan.getReadIndex()], cigarBuffer_);
        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, orphan.highClipped, shadowList[orphan.getReadIndex()], cigarBuffer_);
        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, shadow.highClipped, shadowList[shadow.getReadIndex()], cigarBuffer_);
        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, shadow.highClipped, orphanList[shadow.getReadIndex()], cigarBuffer_);

//        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, shadow.getReadIndex() ? orphan.highClipped : shadow.highClipped, r0Orphan, cigarBuffer_);
//        peAdapterTrimmer_.trimPEAdapterCycles(contigList, readMetadataList, orphan.getReadIndex() ? orphan.highClipped : shadow.highClipped, r1Orphan, cigarBuffer_);

        // recompute orphan scores after adapter trimming
        // TODO: Note we're not accounting for trimmed rescued shadows of another orphan which is probably not good
        orphan.alignmentScore = computeAlignmentScore(orphan, rog, orphanList[orphan.getReadIndex()]);
        orphan.mapQ = alignmentScoreToMapq(orphan.alignmentScore);

        //repeat all pair building on trimmed fragments, but don't base best pair selection on the results
        //TODO: this is somewhat inefficient as we don't need to accumulate resulting pairs, only the pair probabilities
        buildRescuedPairs(anomalousPairHandicap_, contigList, rog, shadow.getReadIndex() ? orphan : shadow, shadowList[1], tls, bestRescuedPair);
        buildRescuedPairs(anomalousPairHandicap_, contigList, rog, orphan.getReadIndex() ? orphan : shadow, shadowList[0], tls, bestRescuedPair);
    }
}


void TemplateBuilder::retainBestSplitAlignments(
        const reference::ContigList& contigList,
        const flowcell::ReadMetadata& shadowReadMetadata,
        const bool regularIndelsOnly,
        const FragmentMetadata &alternative,
        FragmentMetadataList &shadowList,
        const FragmentMetadataList::const_iterator semialignedBegin,
        const FragmentMetadataList::iterator semialignedEnd,
        Cigar &cigarBuffer) const
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(alternative.getCluster().getId(), "    alt: " << alternative);
    for (FragmentMetadataList::const_iterator semialigned = semialignedBegin; semialignedEnd != semialigned; ++semialigned)
    {
        if (alternative != *semialigned)
        {
            const std::size_t before = cigarBuffer.size();
            bool keep = false;
            ISAAC_ASSERT_MSG(shadowList.size() < shadowList.capacity() - TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING, "Unexpected number of shadows shadowList.size():" << shadowList.size());
            if (splitReadAligner_.resolveConflict(
                contigList, shadowReadMetadata, regularIndelsOnly, cigarBuffer, shadowList, alternative, *semialigned))
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(alternative.getCluster().getId(), "    new split: " << shadowList.back());
                // don't push into the heap if we'll have to pop it out. This significantly limits the number of alignments
                // we will ever push into the heap. The worst case is that for some miracle every alignment has unique
                // smith waterman score (and even that is not possible). But then there is only under a thousand smith waterman
                // scores there can be on a few hundred bases long read. Then, if for each smith waterman score there is
                // a multitude of log probabilities we are still talking about 10K CIGARS of alignments ever pushed into shadowList.
                if (shadowList.size() < shadowList.capacity() - TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING ||
                    // shadowList.back() is now the new split alignment. *semialignedEnd is the worst of the splits. It is valid
                    // because initially &*semialignedEnd == &shadowList.back()
                    // if the shadowList is already full, don't bother to insert anything new unless new candidate is better
                    // than the worst of the ones we've seen already
                    FragmentMetadata::bestGappedLess(shadowList.back(), *semialignedEnd))
                {
                    keep = true;
                    std::push_heap(semialignedEnd, shadowList.end(), &FragmentMetadata::bestGappedLess);
                    // need to keep some room at the end
                    ISAAC_ASSERT_MSG(shadowList.size() <= shadowList.capacity() - TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING, "Unexpected number of shadows shadowList.size():" << shadowList.size());
                    if (shadowList.size() == shadowList.capacity() - TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING)
                    {
                        std::pop_heap(semialignedEnd, shadowList.end(), &FragmentMetadata::bestGappedLess);
                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(alternative.getCluster().getId(), "    removing bad split: " << shadowList.back());
                        shadowList.pop_back();
                    }
                }
                else
                {
                    // heap is full and this one is worse than the worst. Remember to pop it from the list!
                    shadowList.pop_back();
                }
            }
            if (!keep)
            {
                // don't keep cigars of alignments we don't keep
                cigarBuffer.resize(before);
            }
        }
        else
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(alternative.getCluster().getId(), "    ign: " << alternative);
        }
    }
}


/**
 * \brief the non-template part of doubleCheckImproperPair
 */
templateBuilder::AlignmentType TemplateBuilder::doubleCheckImproperPair(
    const reference::ContigList& contigList,
    const RestOfGenomeCorrection& rog,
    const flowcell::ReadMetadataList& readMetadataList,
    const TemplateLengthStatistics& tls,
    const bool r0Found, const bool r1Found,
    FragmentMetadataLists &fragments,
    FragmentMetadataLists &shadowList,
    BamTemplate& improperTemplate) const
{
    if (r0Found || r1Found)
    {
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadowList_[0].front().getCluster().getId(),"doubleCheckImproperPair: r0Found: " << r0Found);
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadowList_[0].front().getCluster().getId(),"doubleCheckImproperPair: r1Found: " << r1Found);
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadowList_[0].front().getCluster().getId(),"doubleCheckImproperPair: shadowList_[0].front(): " << shadowList_[0].front() << " " << shadowList_[0].front().getStrandReferencePosition());
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadowList_[0].front().getCluster().getId(),"doubleCheckImproperPair: bamTemplate.getFragmentMetadata(0): " << bamTemplate.getFragmentMetadata(0) << " " << bamTemplate.getFragmentMetadata(0).getFStrandReferencePosition());
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadowList_[0].front().getCluster().getId(),"doubleCheckImproperPair: shadowList_[1].front(): " << shadowList_[1].front() << " " << shadowList_[1].front().getStrandReferencePosition());
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(shadowList_[0].front().getCluster().getId(),"doubleCheckImproperPair: bamTemplate.getFragmentMetadata(1): " << bamTemplate.getFragmentMetadata(1) << " " << bamTemplate.getFragmentMetadata(1).getFStrandReferencePosition());

        if (r0Found && r1Found && shadowList[0].front().splitAlignment && shadowList[1].front().splitAlignment &&
            improperTemplate.getFragmentMetadata(0).getStrandReferencePosition() == shadowList[0].front().getStrandReferencePosition() &&
            improperTemplate.getFragmentMetadata(1).getStrandReferencePosition() == shadowList[1].front().getStrandReferencePosition()
            )
        {
            improperTemplate = BamTemplate(shadowList[0].front(), shadowList[1].front(), improperTemplate.isProperPair());
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(improperTemplate.getCluster().getId(),"doubleCheckImproperPair: rescued orphans on both sides. Assuming both mates cross same breakpoint: " << improperTemplate);
        }
        else
        {
            BamTemplate rescuedTemplate(improperTemplate);
            bestRescuedPair_.removeRepeatDuplicates();
            pickRandomRepeatAlignment(improperTemplate.getCluster().getId(), bestRescuedPair_, rescuedTemplate);
            if (!improperTemplate.isBetterThan(anomalousPairHandicap_, rog.getRogCorrection(), rescuedTemplate))
            {
                // figure out which orphan ended up producing best rescued pair
                const unsigned orphanIndex = improperTemplate.getFragmentMetadata(1) == rescuedTemplate.getFragmentMetadata(1);
                ISAAC_ASSERT_MSG(improperTemplate.getFragmentMetadata(orphanIndex) == rescuedTemplate.getFragmentMetadata(orphanIndex),
                                 "Orphan alignments must match between original and rescued:\n" << improperTemplate << "\n" << rescuedTemplate);

                FragmentMetadata orphan = rescuedTemplate.getFragmentMetadata(orphanIndex);
                // score the assumed orphan using both seed and rescued candidates
                scoreAnomalousEnd(rog, fragments[orphanIndex], shadowList[orphanIndex], orphan);

                FragmentMetadata shadow = rescuedTemplate.getFragmentMetadata(!orphanIndex);
                trimDoublecheckedPairPEAdapaters(
                        contigList, rog, readMetadataList, tls,
                        orphan, shadow,
                        improperTemplate.getFragmentMetadata(0), improperTemplate.getFragmentMetadata(1),
                        fragments, shadowList, bestRescuedPair_);

                bestRescuedPair_.appendPairProbability(improperTemplate.getFragmentMetadata(orphanIndex), improperTemplate.getFragmentMetadata(!orphanIndex), false);

                improperTemplate = BamTemplate(orphan, shadow, rescuedTemplate.isProperPair());

                scoreRescuedShadowTemplate(rog, orphan.getReadIndex(), improperTemplate, bestRescuedPair_);

                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(),"doubleCheckImproperPair: rescued  template: " << improperTemplate);
                // Don't fall through. The template scored.
                return templateBuilder::Normal;
            }
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(improperTemplate.getCluster().getId(),"doubleCheckImproperPair: best rescued:\n" << rescuedTemplate <<
                "\nis worse than improper:\n" << improperTemplate << " mq:" << computeAlignmentScore(rog.getRogCorrection(), exp(improperTemplate.getLogProbability()), exp(rescuedTemplate.getLogProbability())));
        }

        FragmentMetadata read0 = improperTemplate.getFragmentMetadata(0);
        FragmentMetadata read1 = improperTemplate.getFragmentMetadata(1);
        // score each end individually
        scoreAnomalousEnd(rog, fragments[0], shadowList[0], read0);
        scoreAnomalousEnd(rog, fragments[1], shadowList[1], read1);
        improperTemplate = BamTemplate(read0, read1, improperTemplate.isProperPair(), improperTemplate.getAlignmentScore());
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(improperTemplate.getCluster().getId(),"doubleCheckImproperPair: recomputed anomalous end scores: " << improperTemplate);
    }
    else if (dodgySingleton(improperTemplate.getFragmentMetadata(0)) || dodgySingleton(improperTemplate.getFragmentMetadata(1)))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(improperTemplate.getCluster().getId(), "doubleCheckImproperPair nothing found, one or both ends are dodgy singletons. Flagging dodgy ");
        return flagDodgyTemplate(improperTemplate);
    }
    else
    {
        // nothing rescued, just stay with suspicious pair.
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(improperTemplate.getCluster().getId(), "doubleCheckImproperPair nothing better found ");
    }

    return templateBuilder::Normal;
}

} // namespace alignment
} // namespace isaac
