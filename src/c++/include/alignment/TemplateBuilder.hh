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
 ** \file TemplateBuilder.hh
 **
 ** \brief Construction of BamTemplate instances
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH

#include <boost/integer/static_min_max.hpp>
#include <boost/noncopyable.hpp>

#include "alignment/BamTemplate.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/Match.hh"
#include "alignment/RestOfGenomeCorrection.hh"
#include "alignment/matchSelector/DebugStorage.hh"
#include "alignment/templateBuilder/FragmentBuilder.hh"
#include "alignment/templateBuilder/ShadowAligner.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "flowcell/ReadMetadata.hh"
#include "reference/Contig.hh"
#include "templateBuilder/BestPairInfo.hh"
#include "templateBuilder/PeAdapterTrimmer.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Utility component creating Template instances from Seed Matches.
 **
 ** The intended use is to create an instance of a TemplateBuilder for each
 ** thread and delegate to that instance the identification of the most likely
 ** template for each cluster. This is done by invoking the buildTemplate method.
 **/
class TemplateBuilder: boost::noncopyable
{
    static const unsigned READS_IN_A_PAIR = 2;
public:
    typedef std::array<FragmentMetadataList, READS_IN_A_PAIR> FragmentMetadataLists;
    typedef short DodgyAlignmentScore;
    static const DodgyAlignmentScore DODGY_ALIGNMENT_SCORE_UNKNOWN=255;
    static const DodgyAlignmentScore DODGY_ALIGNMENT_SCORE_UNALIGNED=-1;
    static const unsigned POSSIBLY_SEMIALIGNED_MISMATCH_COUNT = 8;
    //static PairClassCounts pairClassCounts_;

    TemplateBuilder(
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
        const unsigned anomalousPairScoreMin,
        const bool reserveBuffers);

    const FragmentMetadataLists &getFragments() const {return candidates_;}

    template <typename MatchFinderT>
    templateBuilder::AlignmentType buildTemplate(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList &readMetadataList,
        const SequencingAdapterList &sequencingAdapters,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const bool withGaps,
        const MatchFinderT &matchFinder,
        BamTemplate &bamTemplate);

    /**
     * \brief public for TemplateDetector
     */
    template <typename MatchFinderT>
    templateBuilder::AlignmentType buildFragments(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const unsigned seedRepeatThreshold,
        const Cluster &cluster,
        const bool withGaps);

    /**
     * \brief public for Unit tests
     */
    templateBuilder::AlignmentType buildCombinationTemplate(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const flowcell::ReadMetadataList &readMetadataList,
        FragmentMetadataLists &fragments,
        const Cluster &cluster,
        const TemplateLengthStatistics &tls,
        BamTemplate &bamTemplate) const;

private:
    // BEST_SHADOWS_TO_KEEP includes semialigned shadows and SV candidates.
    // higher numbers lead to too many candidates stored causing cigar buffer running out of capacity
    static const std::size_t BEST_SHADOWS_TO_KEEP = 1000;
    static const std::size_t TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING = 10;
    static const unsigned SEMIALIGNED_MATCHES_MIN = 16;

    const unsigned repeatThreshold_;
    const unsigned seedLength_;
    const unsigned matchFinderTooManyRepeats_;
    const unsigned matchFinderWayTooManyRepeats_;
    const unsigned matchFinderShadowSplitRepeats_;

    const bool scatterRepeats_;
    const bool rescueShadows_;
    const bool anchorMate_;
    const DodgyAlignmentScore dodgyAlignmentScore_;
    const unsigned anomalousPairScoreMin_;

    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const unsigned smitWatermanGapsMax_;
    const bool splitAlignments_;
    const AlignmentCfg &alignmentCfg_;
    templateBuilder::PeAdapterTrimmer peAdapterTrimmer_;

    /// Buffer for the cigar strings of aligned and rescued reads
    mutable Cigar cigarBuffer_;
    /**
     ** \brief All FragmentMetadata for all reads
     **
     ** candidates_[i] is the list of fragments for read i.
     **/
    FragmentMetadataLists candidates_;
    /// Helper component to align fragments individually
    const templateBuilder::FragmentBuilder fragmentBuilder_;

    static const unsigned SHADOW_ALIGNER_KMER_LENGTH = 8;
    mutable templateBuilder::ShadowAligner<SHADOW_ALIGNER_KMER_LENGTH> shadowAligner_;
    mutable templateBuilder::SplitReadAligner splitReadAligner_;

    /// Buffer for the list of shadows rescued by the shadow aligner
    mutable FragmentMetadataLists shadowList_;

    /// Holds the information about pairs obtained via combining alignments from MatchFinder
    mutable templateBuilder::BestPairInfo bestCombinationPairInfo_;
    /// Holds the information about the pairs rescued via rescueShadow or buildDisjoinedTemplate
    mutable templateBuilder::BestPairInfo bestRescuedPair_;

    template <typename MatchFinderT>
    templateBuilder::AlignmentType buildTemplateFromSeeds(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList &readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const bool withGaps,
        const MatchFinderT &matchFinder,
        BamTemplate &bamTemplate);

    template <typename MatchFinderT>
    templateBuilder::AlignmentType rescueSingletonShadow(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList &readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const MatchFinderT &matchFinder,
        FragmentMetadataLists &fragments,
        BamTemplate &bamTemplate) const;

    template <typename MatchFinderT>
    void rescueBestOrphansShadowTemplates(
        const reference::ContigList& contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList& readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
        const MatchFinderT &matchFinder,
        const TemplateLengthStatistics& tls,
        const FragmentMetadataList& orphans,
        FragmentMetadataList& shadowList,
        templateBuilder::BestPairInfo &bestRescuedPair) const;

            template <typename MatchFinderT>
    templateBuilder::AlignmentType rescueAnomalousRepeatMate(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList &readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const MatchFinderT &matchFinder,
        FragmentMetadataLists &fragments,
        BamTemplate &bamTemplate);

    template <typename MatchFinderT>
    void fixSemialignedMates(
            const reference::ContigList &contigList,
            const RestOfGenomeCorrection &rog,
            const flowcell::ReadMetadataList &readMetadataList,
            templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
            const Cluster &cluster,
            const TemplateLengthStatistics &templateLengthStatistics,
            const MatchFinderT &matchFinder,
            FragmentMetadataLists &fragments,
            BamTemplate &bamTemplate);

    /// Helper method to select the best fragment for single-ended runs
    void pickBestFragment(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const FragmentMetadataList &fragmentList,
        BamTemplate &result) const;

    bool pickBestPair(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const RestOfGenomeCorrection &rog,
        FragmentMetadataLists &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        BamTemplate &bamTemplate) const;

    bool locateBestAnchoredPair(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &rog,
        const FragmentMetadataLists &fragments,
        const TemplateLengthStatistics &templateLengthStatistics,
        unsigned bestPairsMax,
        templateBuilder::BestPairInfo &ret) const;

    void pickRandomRepeatAlignment(
        const unsigned clusterId,
        const templateBuilder::BestPairInfo &bestPair,
        BamTemplate &bamTemplate) const;

    void scoreBestPair(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const FragmentMetadataLists &fragments,
        templateBuilder::BestPairInfo &bestPair,
        BamTemplate &bamTemplate) const;

    void buildSingletonShadowTemplate(
        const RestOfGenomeCorrection &restOfGenomeCorrection,
        const TemplateLengthStatistics &templateLengthStatistics,
        const FragmentMetadataList  &orphans,
        const FragmentMetadata &bestOrphan,
        const unsigned orphanIndex,
        BamTemplate &bamTemplate) const;

    /// Helper method to select the fragment with the highest logProbability
    FragmentMetadataList::const_iterator getBestFragment(
        const FragmentMetadataList &fragmentList) const;

    FragmentMetadataList::iterator pushSemialignedDown(
        const reference::ContigList& contigList,
        FragmentMetadataList& shadowList,
        std::size_t &semialignedHeadLength) const;

    void retainBestSplitAlignments(
            const reference::ContigList& contigList,
            const flowcell::ReadMetadata& shadowReadMetadata,
            const bool regularIndelsOnly,
            const FragmentMetadata &alternative,
            FragmentMetadataList &shadowList,
            const FragmentMetadataList::const_iterator semialignedBegin,
            const FragmentMetadataList::iterator semialignedEnd,
            Cigar &cigarBuffer) const;

    template <typename MatchFinderT>
    bool searchForStructuralVariant(
        const reference::ContigList& contigList,
        const flowcell::ReadMetadata& shadowReadMetadata,
        const bool regularIndelsOnly,
        const unsigned filterContigId,
        templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
        const MatchFinderT &matchFinder,
        const Cluster& cluster,
        FragmentMetadataList &shadowList,
        Cigar &cigarBuffer) const;

    void trimShadowPairPEAdapaters(
        const reference::ContigList& contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList& readMetadataList,
        const TemplateLengthStatistics& tls,
        FragmentMetadata &orphan,
        FragmentMetadata &shadow,
        FragmentMetadataList &orphanList,
        FragmentMetadataList &shadowList,
        templateBuilder::BestPairInfo &bestRescuedPair) const;

    void trimDoublecheckedPairPEAdapaters(
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
        templateBuilder::BestPairInfo &bestRescuedPair) const;


    static void scoreAnomalousEnd(
        const RestOfGenomeCorrection &rog,
        const FragmentMetadataList &seedCandidates,
        FragmentMetadataList &shadowCandidates,
        FragmentMetadata &fragment);

    template <typename MatchFinderT>
    bool rescueMate(
        const reference::ContigList& contigList,
        const bool sameContigOnly,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList& readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const unsigned orphanIndex,
        const TemplateLengthStatistics& templateLengthStatistics,
        FragmentMetadataList& orphanCandidates,
        FragmentMetadataList& shadowList,
        BamTemplate &bamTemplate) const;

    template <typename MatchFinderT>
    bool rescueShadowTemplate(
        const reference::ContigList& contigList,
        const bool sameContigOnly,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList& readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const unsigned orphanIndex,
        const TemplateLengthStatistics& templateLengthStatistics,
        FragmentMetadataList& orphanCandidates,
        FragmentMetadataList& shadowList,
        BamTemplate &bamTemplate) const;

    template <typename MatchFinderT>
    bool rescueShadows(
        const reference::ContigList& contigList,
        const bool sameContigOnly,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList& readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const FragmentMetadata& orphan,
        const TemplateLengthStatistics& tls,
        FragmentMetadataList& shadowList,
        templateBuilder::BestPairInfo &bestRescuedPair) const;

    templateBuilder::AlignmentType doubleCheckImproperPair(
        const reference::ContigList& contigList,
        const RestOfGenomeCorrection& rog,
        const flowcell::ReadMetadataList& readMetadataList,
        const TemplateLengthStatistics& tls,
        const bool r0Found, const bool r1Found,
        FragmentMetadataLists &fragments,
        FragmentMetadataLists &shadowList,
        BamTemplate& bamTemplate) const;

    template <typename MatchFinderT>
    templateBuilder::AlignmentType doubleCheckImproperPair(
        const reference::ContigList& contigList,
        const RestOfGenomeCorrection& rog,
        const flowcell::ReadMetadataList& readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const TemplateLengthStatistics& tls,
        FragmentMetadataLists &fragments,
        BamTemplate& improperTemplate);

    templateBuilder::AlignmentType flagDodgyTemplate(BamTemplate &bamTemplate) const;
};


inline unsigned computeAlignmentScore(
    const FragmentMetadata &fragment,
    const RestOfGenomeCorrection &restOfGenomeCorrection,
    const FragmentMetadataList &fragmentList)
{
//    ISAAC_ASSERT_MSG(fragmentList.end() == std::adjacent_find(fragmentList.begin(), fragmentList.end()),
//                     "Expecting unique alignments in fragment list got\n" << *std::adjacent_find(fragmentList.begin(), fragmentList.end()) <<
//                     " and\n" << *(std::adjacent_find(fragmentList.begin(), fragmentList.end())+1));

    bool bestIgnored = false;
    double neighborProbability = 0.0;
    for (const FragmentMetadata &f : fragmentList)
    {
        if (!bestIgnored && f.logProbability == fragment.logProbability)
        {
            bestIgnored = true;
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore neighborProbability:" << neighborProbability << " skip " << f);
        }
        else
        {
            neighborProbability += exp(f.logProbability);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore neighborProbability:" << neighborProbability << "  add " << f);
        }
    }

    const unsigned ret = computeAlignmentScore(
        restOfGenomeCorrection.getReadRogCorrection(fragment.getReadIndex()),
        exp(fragment.logProbability),
        neighborProbability);

// -ffast-math produces 31 with floor(32)...
//    const double rog = restOfGenomeCorrection.getReadRogCorrection(fragment.getReadIndex());
//    const double alignmentProbability = exp(fragment.logProbability);
//    const double otherAlignmentsProbability = neighborProbability;
//
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore exp(fragment.logProbability):" << alignmentProbability);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore restOfGenomeCorrection.getReadRogCorrection(fragment.getReadIndex()):" << rog);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore neighborProbability:" << otherAlignmentsProbability);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore (otherAlignmentsProbability + alignmentProbability + restOfGenomeCorrection):" << (otherAlignmentsProbability + alignmentProbability + rog));
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore (otherAlignmentsProbability + rog):" << (otherAlignmentsProbability + rog));
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore div:" << (otherAlignmentsProbability + rog) / (otherAlignmentsProbability + alignmentProbability + rog));
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore log10(div):" << log10((otherAlignmentsProbability + rog) / (otherAlignmentsProbability + alignmentProbability + rog)));
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore -10.0 * log10(div):" << -10.0*log10((otherAlignmentsProbability + rog) / (otherAlignmentsProbability + alignmentProbability + rog)));
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragment.getCluster().getId(), "computeAlignmentScore floor(-10.0 * log10(div)):" << floor(-10.0*log10((otherAlignmentsProbability + rog) / (otherAlignmentsProbability + alignmentProbability + rog))));
    return ret;
}

template <typename MatchFinderT>
bool TemplateBuilder::searchForStructuralVariant(
    const reference::ContigList& contigList,
    const flowcell::ReadMetadata& shadowReadMetadata,
    const bool regularIndelsOnly,
    const unsigned filterContigId,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    const MatchFinderT &matchFinder,
    const Cluster& cluster,
    FragmentMetadataList &shadowList,
    Cigar &cigarBuffer) const
{
    const bool bestWasGapped = shadowList.front().gapCount;
    const FragmentMetadataList::iterator firstSplit = shadowList.end();

    if (shadowList.size() >= shadowList.capacity() - TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "   no capacity to split shadows");
        return false;
    }

    std::size_t semialignedHeadLength = 0;
    const FragmentMetadataList::const_iterator firstSemialigned =
        pushSemialignedDown(contigList, shadowList, semialignedHeadLength);
    if (shadowList.end() == firstSemialigned)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "   no semialigned shadows");
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
            cluster.getId(), "   found " << std::distance<FragmentMetadataList::const_iterator>(firstSemialigned, shadowList.end()) << " semialigned shadows, semialignedHeadLength:" << semialignedHeadLength);

        if (semialignedHeadLength > seedLength_ &&
            templateBuilder::Normal != fragmentBuilder_.buildAllHeadAnchored(
            contigList, shadowReadMetadata, filterContigId, matchFinderShadowSplitRepeats_,
            // the first offset from which a seed cannot be made
            std::min<unsigned>(shadowReadMetadata.getLength() - seedLength_, semialignedHeadLength),
            adapterClipper, matchFinder, cluster,
            [this, &contigList, &firstSemialigned, &firstSplit, &shadowReadMetadata, regularIndelsOnly, &cigarBuffer, &shadowList](
                    FragmentMetadata &&alternative)
            {
                if (!clippedByReference(contigList, alternative) && alternative.recomputeHeadAnchor(contigList))
                {
                    retainBestSplitAlignments(
                            contigList, shadowReadMetadata, regularIndelsOnly, alternative,
                            shadowList, firstSemialigned, firstSplit, cigarBuffer);
                }
                else
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(alternative.getCluster().getId(), "    no head anchor: " << alternative);
                }
            }))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "    could not build any alignments from seeds with filterContigId:" << filterContigId);
        }
    }

    if (shadowList.end() != firstSplit)
    {
        return putBestOnTop<true>(shadowList);
    }
    // this needs to be done because even if searchForStructuralVariant does not find anything, it rearranges
    // shadow list so that the semialigned candidates are at the bottom.
    else if (bestWasGapped)
    {
        putBestOnTop<true>(shadowList);
    }
    else
    {
        putBestOnTop<false>(shadowList);
    }
    return false;
}

inline void upateBestRescuedPair(
    const unsigned anomalousPairScoreMin,
    const RestOfGenomeCorrection &rog,
    const TemplateLengthStatistics& tls,
    const FragmentMetadata& orphan,
    const FragmentMetadata& rescuedShadow,
    templateBuilder::BestPairInfo& ret)
{
    const isaac::alignment::TemplateLengthStatistics::CheckModelResult model = tls.checkModel(orphan, rescuedShadow);
    const bool properPair = TemplateLengthStatistics::Nominal == model || TemplateLengthStatistics::Undersized == model;
    const templateBuilder::PairInfo pairInfo(orphan, rescuedShadow, properPair);

    // Notice that all pairs we deal with here are properly oriented as this is how the rescue works. Some of them are
    // oversized, but we allow oversized pairs to be added to the proper pair probabilities list since this is the
    // reason we rescue outside insert size. However we don't allow them to be marked as proper. This way they reduce
    // proper pair alignment score, but get scored with the chimera and splits if they happen to be the best choice.
    const bool hasSplits = orphan.isSplit() || rescuedShadow.isSplit();

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "upateBestRescuedPair:\n" <<
        "orphan:" << orphan << "\n" <<
        "shadow:" << rescuedShadow << "\n" <<
        "ret:" << ret);

//    if (!ret..isBetterThan(anomalousPairScoreMin_, rog.getRogCorrection(), rescuedTemplate))
    if (ret.isWorseThan(anomalousPairScoreMin, rog, pairInfo))
    {
        ret.resetBest(pairInfo, orphan, rescuedShadow, !hasSplits, properPair);
    }
    else if (ret.isAsGood(pairInfo))
    {
        ret.appendBest(orphan, rescuedShadow, !hasSplits, properPair);
    }
    else
    {
        ret.appendPairProbability(orphan, rescuedShadow, !hasSplits);
    }
}

inline int updateBestAnchoredPair(
    const unsigned anomalousPairScoreMin,
    const RestOfGenomeCorrection &rog,
    const TemplateLengthStatistics &tls,
    const FragmentMetadata &fragmentA,
    const FragmentMetadata &fragmentB,
    templateBuilder::BestPairInfo &ret)
{
    const isaac::alignment::TemplateLengthStatistics::CheckModelResult model = tls.checkModel(fragmentA, fragmentB);
    const bool properPair = TemplateLengthStatistics::Nominal == model || TemplateLengthStatistics::Undersized == model;
    const templateBuilder::PairInfo pairInfo(fragmentA, fragmentB, properPair);

    // unlike rescued pairs, there is no extra allowance for oversized pairs when global alignment is performed.
    // The logic is that cases resolved by seeds are fairly trivial, so it is unlikely that multiple alignments
    // within same location are similar enough for insert size to be a criterion worth coding for.
    if (ret.isWorseThan(anomalousPairScoreMin, rog, pairInfo))
    {
        ret.resetBest(pairInfo, fragmentA, fragmentB, properPair, properPair);
        return -1;
    }
    else if (ret.isAsGood(pairInfo))
    {
        ret.appendBest(fragmentA, fragmentB, properPair, properPair);
        return 1;
    }
    else
    {
        ret.appendPairProbability(fragmentA, fragmentB, properPair);
    }
    return 0;
}

inline void buildRescuedPairs(
    const unsigned anomalousPairScoreMin,
    const reference::ContigList& contigList,
    const RestOfGenomeCorrection &rog,
    const FragmentMetadata& orphan,
    const FragmentMetadataList &shadowList,
    const TemplateLengthStatistics& tls,
    templateBuilder::BestPairInfo &bestRescuedPair)
{
    for(const FragmentMetadata &rescuedShadow : shadowList)
    {
        //all rescued pairs match model and are proper pairs even if they have splits because that's how we found them
        upateBestRescuedPair(anomalousPairScoreMin, rog, tls, orphan, rescuedShadow, bestRescuedPair);
        bestRescuedPair.appendSingleProbability(rescuedShadow);
    }
}

inline bool isVeryBadShadow(
    const double rog,
    const FragmentMetadata &bestRescuedShadow)
{
    static const unsigned BAD_SHADOW_NO_NEIGHBOR_ALIGNMENT_SCORE = 7;
    const unsigned noNeighborAlignmentScore = computeAlignmentScore(rog, exp(bestRescuedShadow.logProbability), 0.0);
    if (BAD_SHADOW_NO_NEIGHBOR_ALIGNMENT_SCORE > noNeighborAlignmentScore)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bestRescuedShadow.getCluster().getId(), "Could not improve a semialigned shadow " << bestRescuedShadow << " bad standalone mapq:" << noNeighborAlignmentScore);
        return true;
    }

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bestRescuedShadow.getCluster().getId(), "Keeping semialigned shadow " << bestRescuedShadow << " exp(bestRescuedShadow.logProbability):" << exp(bestRescuedShadow.logProbability));
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bestRescuedShadow.getCluster().getId(), "Keeping semialigned shadow " << bestRescuedShadow << " rog:" << rog);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bestRescuedShadow.getCluster().getId(), "Keeping semialigned shadow " << bestRescuedShadow << " standalone mapq:" << noNeighborAlignmentScore);
    return false;
}


inline void TemplateBuilder::scoreAnomalousEnd(
    const RestOfGenomeCorrection &rog,
    const FragmentMetadataList &seedCandidates,
    FragmentMetadataList &shadowCandidates,
    FragmentMetadata &fragment)
{
    // mix everything in
    shadowCandidates.insert(shadowCandidates.end(), seedCandidates.begin(), seedCandidates.begin() + std::min(seedCandidates.size(), TOP_BEST_SEED_CANDIDATES_FOR_ANOMALOUS_SCORING));
    // sort to remove any occasional duplicate alignments
    std::sort(shadowCandidates.begin(), shadowCandidates.end());
    shadowCandidates.erase(std::unique(shadowCandidates.begin(), shadowCandidates.end()), shadowCandidates.end());

    const unsigned newScore = computeAlignmentScore(fragment, rog, shadowCandidates);
    if (newScore < fragment.alignmentScore)
    {
        fragment.alignmentScore = newScore;
        fragment.mapQ = alignmentScoreToMapq(fragment.alignmentScore);
    }
}

inline void scoreRescuedShadowTemplate(
    const RestOfGenomeCorrection &rog,
    const unsigned orphanIndex,
    BamTemplate &bamTemplate,
    templateBuilder::BestPairInfo &bestRescuedPair)
{
    FragmentMetadata &orphan = bamTemplate.getFragmentMetadata(orphanIndex);
    FragmentMetadata &shadow = bamTemplate.getFragmentMetadata(!orphanIndex);

    const double otherTemplateProbability =
        bestRescuedPair.sumUniquePairProbabilities(
            orphan.logProbability + shadow.logProbability,
            bestRescuedPair.repeatsCount(),
            bamTemplate.isProperPair());
    bamTemplate.setAlignmentScore(computeAlignmentScore(rog.getRogCorrection(), bestRescuedPair.probability(), otherTemplateProbability));
    ISAAC_ASSERT_MSG(UNKNOWN_MAPQ != orphan.mapQ, "Invalid orphan.alignmentScore: " << orphan);
    shadow.mapQ = pickMapQFromMate(orphan.mapQ, bamTemplate.getAlignmentScore());
}


inline bool dodgySingleton(const FragmentMetadata &singleton)
{
    return singleton.uncheckedSeeds &&
        (singleton.gapCount ||
            // this can be improved. For the moment the most basic approach is that if a repeat-seed singleton has 1 mismatch,
            // the only locations we could have missed are another bunch of 1-mismatch alignments with mismatches in the seeds we've visited
            // 1 < singleton.mismatchCount

            // Above is a nice idea but then it does lead to high-confidence misplaced pairs when true alignment has 1 or more mismatches.
            // Considering any imperfect singleton to be dodgy for now.
            singleton.mismatchCount);
}

/**
 ** \brief Attempt to rescue shadow of a single orphan
 **
 **/
template <typename MatchFinderT>
templateBuilder::AlignmentType TemplateBuilder::rescueSingletonShadow(
    const reference::ContigList &contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList &readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics,
    const MatchFinderT &matchFinder,
    FragmentMetadataLists &fragments,
    BamTemplate &bamTemplate) const
{
    const unsigned orphanIndex = fragments[0].empty();
    const FragmentMetadata &orphan = bamTemplate.getFragmentMetadata(orphanIndex);

    if (orphan.decoyAlignment)
    {
        return flagDodgyTemplate(bamTemplate);
    }
    else if (orphan.isUniquelyAligned() && !dodgySingleton(orphan))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Unique orphan: " << orphan << " Rescuing shadow...");
        rescueShadowTemplate(
            contigList, true, rog, readMetadataList, adapterClipper, matchFinder,
            orphanIndex, templateLengthStatistics, fragments[orphanIndex], shadowList_[!orphanIndex], bamTemplate);
        // regardless of success return Normal as in worst case we have singleton/shadow template
    }

    return templateBuilder::Normal;
}

/**
 * \brief run through best orphans and rescue their shadows as well. This is meant to be called for repeat mate
 *        to make sure that we don't score unique mate/rescued shadow pair too high
 */
template <typename MatchFinderT>
void TemplateBuilder::rescueBestOrphansShadowTemplates(
    const reference::ContigList& contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList& readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    const MatchFinderT &matchFinder,
    const TemplateLengthStatistics& tls,
    const FragmentMetadataList& orphans,
    FragmentMetadataList& shadowList,
    templateBuilder::BestPairInfo &bestRescuedPair) const
{
    ISAAC_ASSERT_MSG(!bestRescuedPair_.empty(), "Expected to have a rescued pair of the unique mate");
    const FragmentMetadata &bestOrphan = orphans.front();
    for (const FragmentMetadata &repeatMate : orphans)
    {
        if (!bestOrphan.isBetterUngapped(repeatMate))
        {
            const std::size_t before = cigarBuffer_.size();
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeatMate.getCluster().getId(), "rescueBestOrphansShadowTemplates for repeatMate: " << repeatMate);
            // we're not interested in shadows, only in updating pair probability
            rescueShadows(contigList, true, rog, readMetadataList, adapterClipper,
                matchFinder, repeatMate, tls, shadowList, bestRescuedPair);
            // cigar buffer is not meant to hold all those shadow cigars. And we don't need them.
            cigarBuffer_.resize(before);
        }
    }
}

/**
 ** \brief Attempt to rescue shadow of a single orphan
 **
 **/
template <typename MatchFinderT>
templateBuilder::AlignmentType TemplateBuilder::rescueAnomalousRepeatMate(
    const reference::ContigList &contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList &readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const Cluster &cluster,
    const TemplateLengthStatistics &tls,
    const MatchFinderT &matchFinder,
    FragmentMetadataLists &fragments,
    BamTemplate &bamTemplate)
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Pair with low-scored mate. Checking that pair is unique within rescue range: " << bamTemplate);
    // to rescue the low-scored mate we need to make sure the pair is unique within the template rescue range
    const unsigned uniqueMateIndex = !bamTemplate.getFragmentMetadata(0).isUniquelyAligned();
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "One end is not unique. Rescuing: " << bamTemplate);
    BamTemplate uniqueOrphanRescued(bamTemplate);
    FragmentMetadata &orphan = bamTemplate.getFragmentMetadata(uniqueMateIndex);
    bestRescuedPair_.clear();
    if (rescueShadows(
        contigList,
        true,
        rog,
        readMetadataList,
        adapterClipper,
        matchFinder,
        orphan,
        tls,
        shadowList_[!uniqueMateIndex],
        bestRescuedPair_))
    {
        bestRescuedPair_.removeRepeatDuplicates();
        pickRandomRepeatAlignment(orphan.getCluster().getId(), bestRescuedPair_, uniqueOrphanRescued);
        if (!bamTemplate.isBetterThan(anomalousPairScoreMin_, rog.getRogCorrection(), uniqueOrphanRescued))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(), "Rescued template:\n" << uniqueOrphanRescued <<
                "is same or better:\n" << bamTemplate <<
                " : " << computeAlignmentScore(rog.getRogCorrection(), exp(bamTemplate.getLogProbability()), exp(uniqueOrphanRescued.getLogProbability())) <<
                " >= " << anomalousPairScoreMin_);

            // check if rescuing from the other end ruins our bubble
            rescueBestOrphansShadowTemplates(
                contigList,
                rog, readMetadataList, adapterClipper, matchFinder,
                tls, fragments[!uniqueMateIndex], shadowList_[uniqueMateIndex], bestRescuedPair_);

            // trim, score and stay positive
            FragmentMetadata &shadow = uniqueOrphanRescued.getFragmentMetadata(!uniqueMateIndex);
            trimShadowPairPEAdapaters(
                contigList, rog, readMetadataList, tls, orphan, shadow, fragments[uniqueMateIndex], shadowList_[!uniqueMateIndex], bestRescuedPair_);

            scoreRescuedShadowTemplate(rog, uniqueMateIndex, uniqueOrphanRescued, bestRescuedPair_);

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(), "Rescued template:\n" << uniqueOrphanRescued <<
                "after repeat mate rescue:\n" << bamTemplate <<
                " : " << computeAlignmentScore(rog.getRogCorrection(), exp(bamTemplate.getLogProbability()), exp(uniqueOrphanRescued.getLogProbability())) <<
                " >= " << anomalousPairScoreMin_);

            bamTemplate = uniqueOrphanRescued;
        }
        else if (dodgySingleton(bamTemplate.getFragmentMetadata(uniqueMateIndex)))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(),
                                                   "nothing good found when trying to rescue repeat end of an improper pair with dodgy singleton. Flagging dodgy ");
            return flagDodgyTemplate(bamTemplate);
        }
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(), "Keeping the seed discovered template: " << bamTemplate);
    }

    return templateBuilder::Normal;
}

/**
 ** \brief Attempt to introduce splits into one or both mates to improve alignment
 */
template <typename MatchFinderT>
void TemplateBuilder::fixSemialignedMates(
        const reference::ContigList &contigList,
        const RestOfGenomeCorrection &rog,
        const flowcell::ReadMetadataList &readMetadataList,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
	    const Cluster &cluster,
        const TemplateLengthStatistics &templateLengthStatistics,
        const MatchFinderT &matchFinder,
        FragmentMetadataLists &fragments,
        BamTemplate &bamTemplate)
{
    // figure out which of the mates is possibly semialigned
    const std::size_t one = bamTemplate.getFragmentMetadata(1).possiblySemialigned(contigList, seedLength_);
    const std::size_t zero = !one;

    // if one that we've already checked or the other one
    if (one || bamTemplate.getFragmentMetadata(one).possiblySemialigned(contigList, seedLength_))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Possibly semialigned r" << one << ". Rescuing : " << bamTemplate);
        BamTemplate rescuedTemplate(bamTemplate);
        if (rescueMate(
                contigList, !splitAlignments_, rog, readMetadataList, adapterClipper, matchFinder,
                zero, templateLengthStatistics, fragments[zero], shadowList_[one], rescuedTemplate) &&
                !bamTemplate.isBetterThan(anomalousPairScoreMin_, rog.getRogCorrection(), rescuedTemplate))
        {
            // add the original pair in so that rescued mate scoring accounts for the original
            bestRescuedPair_.appendPairProbability(bamTemplate.getFragmentMetadata(zero), bamTemplate.getFragmentMetadata(one), true);
            scoreRescuedShadowTemplate(rog, zero, rescuedTemplate, bestRescuedPair_);

            bamTemplate = rescuedTemplate;
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(), "Rescued template is better: " << bamTemplate);
            // if the other mate is possibly semialigned
            if (bamTemplate.getFragmentMetadata(zero).possiblySemialigned(contigList, seedLength_))
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Possibly semialigned r" << zero << ". Rescuing : " << bamTemplate);
                BamTemplate rescuedTemplate(bamTemplate);
                if (rescueMate(
                        contigList, !splitAlignments_, rog, readMetadataList, adapterClipper, matchFinder,
                        one, templateLengthStatistics, fragments[one], shadowList_[zero], rescuedTemplate) &&
                        !bamTemplate.isBetterThan(anomalousPairScoreMin_, rog.getRogCorrection(), rescuedTemplate))
                {
                    bamTemplate = rescuedTemplate;
                    // both ends had to be rescued. Just score them independently usign seed and shadow rescue candidates
                    scoreAnomalousEnd(rog, fragments[zero], shadowList_[zero], bamTemplate.getFragmentMetadata(zero));
                    scoreAnomalousEnd(rog, fragments[one], shadowList_[one], bamTemplate.getFragmentMetadata(one));
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(), "Rescued template is better: " << bamTemplate);
                }
            }
        }
    }
}

template <typename MatchFinderT>
bool TemplateBuilder::rescueMate(
    const reference::ContigList& contigList,
    const bool sameContigOnly,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList& readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    const MatchFinderT &matchFinder,
    const unsigned orphanIndex,
    const TemplateLengthStatistics& tls,
    FragmentMetadataList& orphanCandidates,
    FragmentMetadataList& shadowList,
    BamTemplate &bamTemplate) const
{
    FragmentMetadata &orphan = bamTemplate.getFragmentMetadata(orphanIndex);
    bestRescuedPair_.clear();
    if (rescueShadows(
        contigList,
        sameContigOnly,
        rog,
        readMetadataList,
        adapterClipper,
        matchFinder,
        orphan,
        tls,
        shadowList,
        bestRescuedPair_))
    {
        bestRescuedPair_.removeRepeatDuplicates();
        pickRandomRepeatAlignment(orphan.getCluster().getId(), bestRescuedPair_, bamTemplate);

        FragmentMetadata &shadow = bamTemplate.getFragmentMetadata(!orphan.getReadIndex());
        trimShadowPairPEAdapaters(
            contigList, rog, readMetadataList, tls, orphan, shadow, orphanCandidates, shadowList, bestRescuedPair_);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(),"rescueMate: rescued: " << bamTemplate);

        return true;
    }
    return false;
}

template <typename MatchFinderT>
bool TemplateBuilder::rescueShadowTemplate(
    const reference::ContigList& contigList,
    const bool sameContigOnly,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList& readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    const MatchFinderT &matchFinder,
    const unsigned orphanIndex,
    const TemplateLengthStatistics& tls,
    FragmentMetadataList& orphanCandidates,
    FragmentMetadataList& shadowList,
    BamTemplate &bamTemplate) const
{
    if (rescueMate(
        contigList, sameContigOnly, rog, readMetadataList, adapterClipper,
        matchFinder, orphanIndex, tls, orphanCandidates, shadowList, bamTemplate))
    {
        scoreRescuedShadowTemplate(rog, orphanIndex, bamTemplate, bestRescuedPair_);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(),"rescueShadowTemplate: rescued: " << bamTemplate);
        return true;
    }
    return false;
}


template <typename MatchFinderT>
bool TemplateBuilder::rescueShadows(
    const reference::ContigList& contigList,
    const bool sameContigOnly,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList& readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    const MatchFinderT &matchFinder,
    const FragmentMetadata& orphan,
    const TemplateLengthStatistics& tls,
    FragmentMetadataList& shadowList,
    templateBuilder::BestPairInfo &bestRescuedPair) const
{
    const unsigned orphanIndex = orphan.getReadIndex();

    shadowList.clear();
    if (shadowAligner_.rescueShadows(
        contigList, orphan,
        // each shadow can produce up to 1 gapped alignment
        shadowList.capacity() / 2, shadowList,
        readMetadataList[!orphanIndex], adapterClipper, tls))
    {
        FragmentMetadata &topShadow = shadowList.front();

        if (smitWatermanGapsMax_ && //BandedSmithWaterman::mismatchesCutoff < bestShadow.mismatchCount)
            (
                // Best candidate has enough mismatches to open a gap.
                topShadow.mismatchCount * alignmentCfg_.normalizedMismatchScore_ >= alignmentCfg_.normalizedGapOpenScore_
//                ||
//                // if we have something that amounts to a long gap, do smith waterman anyway
//                bestShadow.smithWatermanScore >= alignmentCfg_.normalizedMaxGapExtendScore_
            ))
        {
            fragmentBuilder_.realignBadUngappedAlignments(contigList, readMetadataList[!orphanIndex], adapterClipper, shadowList);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Best shadow after smith waterman: " << topShadow);
        }

        // see if best rescued is any good ignoring the alternatives
        // const unsigned assumedlyUniqueMapq = computeAlignmentScore(rog.getReadRogCorrection(shadowReadMetadata.getIndex()), exp(bestRescuedShadow.logProbability), 0.0);
        // if (isVeryBadAlignment(bestRescued))
        if (topShadow.possiblySemialigned(contigList, seedLength_))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Rescued shadow too bad: " << topShadow);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    splitAlignments_: " << splitAlignments_);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    sameContigOnly: " << sameContigOnly);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    orphan.getContigId(): " << orphan.getContigId());
            // search for reasonable indels first. Accept them even if SVs with less mismatches are possbile
            if (!searchForStructuralVariant(
                contigList, readMetadataList[!orphanIndex], true, orphan.getContigId(),
                adapterClipper, matchFinder, orphan.getCluster(), shadowList, cigarBuffer_) &&
                splitAlignments_)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    orphan.getContigId(): " << orphan.getContigId());
                // search for everything if splitting is allowed
                searchForStructuralVariant(
                    contigList, readMetadataList[!orphanIndex],
                    false, sameContigOnly ? orphan.getContigId() : MatchFinderT::NO_CONTIG_FILTER,
                    adapterClipper, matchFinder, orphan.getCluster(), shadowList, cigarBuffer_);
            }
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Best shadow after structural variant: " << topShadow);
        }

        if (isVeryBadShadow(rog.getReadRogCorrection(!orphanIndex), topShadow))
        {
            return false;
        }

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "    Best shadow: " << topShadow);

        buildRescuedPairs(anomalousPairScoreMin_, contigList, rog, orphan, shadowList, tls, bestRescuedPair);
        return true;
    }
    else if (!shadowList.empty())
    {
        // The shadow hits a repetitive region next to one of the orphans
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), (boost::format("Shadow rescue hits a repeat. Orphan: %s") % orphan).str());
    }
    else
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(orphan.getCluster().getId(), "No shadows rescued");
    }
    return false;
}


template <typename MatchFinderT>
templateBuilder::AlignmentType TemplateBuilder::doubleCheckImproperPair(
    const reference::ContigList& contigList,
    const RestOfGenomeCorrection& rog,
    const flowcell::ReadMetadataList& readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const MatchFinderT &matchFinder,
    const TemplateLengthStatistics& tls,
    FragmentMetadataLists &fragments,
    BamTemplate& improperTemplate)
{
//    ISAAC_ASSERT_MSG(bamTemplate.isUniquelyAligned(), "Doublechecking for non-unique template is not allowed " << bamTemplate);
    bestRescuedPair_.clear();
    const bool r0Found = rescueShadows(
        contigList, !splitAlignments_, rog, readMetadataList, adapterClipper, matchFinder,
        improperTemplate.getFragmentMetadata(0), tls, shadowList_[1], bestRescuedPair_);

    const bool r1Found = rescueShadows(
        contigList, !splitAlignments_, rog, readMetadataList, adapterClipper, matchFinder,
        improperTemplate.getFragmentMetadata(1), tls, shadowList_[0], bestRescuedPair_);

    return doubleCheckImproperPair(
        contigList, rog, readMetadataList, tls, r0Found, r1Found, fragments, shadowList_, improperTemplate);
}


/**
 * \brief public for TemplateDetector
 */
template <typename MatchFinderT>
templateBuilder::AlignmentType TemplateBuilder::buildFragments(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadataList &readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const MatchFinderT &matchFinder,
    const unsigned seedRepeatThreshold,
    const Cluster &cluster,
    const bool withGaps)
{
    cigarBuffer_.clear();
    templateBuilder::AlignmentType ret = templateBuilder::Nm;
    for (const flowcell::ReadMetadata &readMetadata : readMetadataList)
    {
        candidates_[readMetadata.getIndex()].clear();

        const templateBuilder::AlignmentType readAlignmentType = fragmentBuilder_.buildBest(
                        contigList, readMetadata,
                        seedRepeatThreshold,
                        adapterClipper,
                        matchFinder, cluster, withGaps,
                        candidates_[readMetadata.getIndex()]);

        ret = combineAlignmentTypes(ret, readAlignmentType);
    }

    return ret;
}

template <typename MatchFinderT>
templateBuilder::AlignmentType TemplateBuilder::buildTemplateFromSeeds(
    const reference::ContigList &contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList &readMetadataList,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics,
    const bool withGaps,
    const MatchFinderT &matchFinder,
    BamTemplate &bamTemplate)
{
    templateBuilder::AlignmentType res = buildFragments(
        contigList, readMetadataList, adapterClipper, matchFinder, matchFinderTooManyRepeats_, cluster, withGaps);

    if (res == templateBuilder::Normal && candidates_[0].empty() != candidates_[1].empty())
    {
        //attempt to use more seeds if we have a singleton or a single-ended alignment with mismatches and
        //unchecked seeds in a non-decoy contig
        const FragmentMetadataList &singletonCandidates = candidates_[candidates_[0].empty()];
        const FragmentMetadata &bestSingleton = *getBestFragment(singletonCandidates);
        if (2 == readMetadataList.size() && bestSingleton.decoyAlignment)
        {
            buildCombinationTemplate(
                contigList, rog, readMetadataList,
                candidates_, cluster, templateLengthStatistics, bamTemplate);
            res =  flagDodgyTemplate(bamTemplate);
        }
        else if (dodgySingleton(bestSingleton))
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Singleton had some repeat seeds, trying with higher threshold : " << candidates_[candidates_[0].empty()].front());
            // in case of singleton discovered while ignoring some seed hits (about 1% for regular human WGS), repeat procedure with higher threshold
            res = buildFragments(
                contigList, readMetadataList, adapterClipper, matchFinder, matchFinderWayTooManyRepeats_, cluster, withGaps);
        }
    }

    if (res == templateBuilder::Normal)
    {
        res = buildCombinationTemplate(
            contigList, rog, readMetadataList,
            candidates_, cluster, templateLengthStatistics, bamTemplate);
    }

    return res;
}

/**
 ** \brief Build the most likely template for a single cluster, given a set of candidate alignments
 **
 ** \return Rm means template ended up not having a single read aligned anywhere.
 **
 ** This method will initialize the internal template of the builder.
 **/
template <typename MatchFinderT>
templateBuilder::AlignmentType TemplateBuilder::buildTemplate(
    const reference::ContigList &contigList,
    const RestOfGenomeCorrection &rog,
    const flowcell::ReadMetadataList &readMetadataList,
    const SequencingAdapterList &sequencingAdapters,
    const Cluster &cluster,
    const TemplateLengthStatistics &templateLengthStatistics,
    const bool withGaps,
    const MatchFinderT &matchFinder,
    BamTemplate &bamTemplate)
{
    templateBuilder::FragmentSequencingAdapterClipper adapterClipper(sequencingAdapters);
    templateBuilder::AlignmentType res = buildTemplateFromSeeds(
        contigList, rog, readMetadataList,
        adapterClipper, cluster, templateLengthStatistics,
        withGaps, matchFinder,
        bamTemplate);

    if (res == templateBuilder::Normal)
    {
        if (2 == readMetadataList.size() && rescueShadows_)
        {
            ISAAC_ASSERT_MSG(!candidates_[0].empty() || !candidates_[1].empty(), "Both paired fragment lists are empty");

            if (candidates_[0].empty() || candidates_[1].empty())
            {
                res = rescueSingletonShadow(
                    contigList, rog, readMetadataList, adapterClipper, cluster, templateLengthStatistics,
                    matchFinder, candidates_, bamTemplate);
            }
            else
            {
                const FragmentMetadata &r0Fragment = bamTemplate.getFragmentMetadata(0);
                const FragmentMetadata &r1Fragment = bamTemplate.getFragmentMetadata(1);
                if (r0Fragment.decoyAlignment || r1Fragment.decoyAlignment)
                {
                    res = flagDodgyTemplate(bamTemplate);
                }
                else
                {
                    if(r0Fragment.isUniquelyAligned() != r1Fragment.isUniquelyAligned())
                    {
                        if (!bamTemplate.isProperPair())
                        {
                            res = rescueAnomalousRepeatMate(
                                contigList, rog, readMetadataList, adapterClipper, cluster,
                                templateLengthStatistics,
                                matchFinder, candidates_, bamTemplate);
                        }
                        //else proper pairs will upgrade repeat mate MAPQ if the pair is unique
                    }
                    else if (1 == r0Fragment.repeatCount && 1 == r1Fragment.repeatCount)
                    {
                        // both unique, need to make sure we don't accept some anomaly without extra checking
                        if (!bamTemplate.isProperPair())
                        {
                            // improper pairs get checked from both ends
                            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "Improper pair. Doublechecking : " << bamTemplate);
                            res = doubleCheckImproperPair(
                                contigList, rog, readMetadataList, adapterClipper, matchFinder,
                                templateLengthStatistics, candidates_, bamTemplate);
                        }
                        else
                        {
                            fixSemialignedMates(
                                contigList, rog, readMetadataList, adapterClipper, cluster,
                                templateLengthStatistics, matchFinder, candidates_, bamTemplate);
                        }
                    }
                }
            }
        }
    }
//    bamTemplate.debugClass_ = 0;//pairClassCounts_.countIn(bamTemplate, rog, templateBuilder::Normal != res);
    return res;
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_HH
