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
 ** \file FragmentBuilder.hh
 **
 ** \brief Utility classes for Fragment building and management for several reads
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH

#include <atomic>

#include "alignment/templateBuilder/GappedAligner.hh"
#include "alignment/templateBuilder/UngappedAligner.hh"
#include "alignment/templateBuilder/SplitReadAligner.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/FragmentMetadata.hh"
#include "alignment/HashMatchFinder.hh"
#include "alignment/Match.hh"
#include "alignment/RestOfGenomeCorrection.hh"
#include "alignment/templateBuilder/FragmentSequencingAdapterClipper.hh"
#include "reference/Contig.hh"
#include "flowcell/ReadMetadata.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

enum AlignmentType
{
    // the order reflects the precedence meaning that when merging results of two read alignment attemps,
    // Normal will override any other value and so on.
    Normal,
    Rm,         // one of the seeds exactly mapped to a high repeat or too many neighbors with the same prefix
    Qc,         // all seeds contain Ns, alignment is not possible
    Nm,       // seeds have no match in the reference.
};


inline AlignmentType combineAlignmentTypes(const AlignmentType left, const AlignmentType right)
{
    return std::min(left, right);
}

/**
 ** \brief Utility component creating and scoring all Fragment instances from a
 ** list Seed Matches for a single Cluster (each Read independently).
 **/
class FragmentBuilder: public boost::noncopyable
{
public:
    FragmentBuilder(
        const bool collectMismatchCycles,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned repeatThreshold,
        const unsigned seedRepeatThreshold,
        const unsigned seedLength,
        const unsigned maxSeedsPerMatch,
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const bool avoidSmithWaterman,
        const unsigned smithWatermanGapSizeMax,
        const bool noSmithWaterman,
        const bool splitAlignments,
        const AlignmentCfg &alignmentCfg,
        Cigar &cigarBuffer,
        const bool reserveBuffers);

    ~FragmentBuilder()
    {
//        if (!countsTraced_)
//        {
//            countsTraced_ = true;
//            for (std::size_t x = 0; x < seedMatchCounts_.size(); ++x)
//            {
//                ISAAC_THREAD_CERR << x << "\t" << seedMatchCounts_[x] << std::endl;
//            }
//            for (std::size_t x = 0; x < matchSeedCounts_.size(); ++x)
//            {
//                std::string line = ":";
//                for (std::size_t y = 0; y < matchSeedCounts_[x].size(); ++y)
//                {
//                    line += "\t";
//                    line += std::to_string(matchSeedCounts_[x][y].first);
//                    line += "/";
//                    line += std::to_string(matchSeedCounts_[x][y].second);
//                }
//                ISAAC_THREAD_CERR << x << ":" << line << std::endl;
//            }
//
//            for (std::size_t x = 0; x < matchSeedCounts_.size(); ++x)
//            {
//                ISAAC_THREAD_CERR << x << ": " << alignedSeedCounts_[x] << std::endl;
//            }
//        }
    }

    template <typename MatchFinderT>
    AlignmentType buildBest(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const std::size_t seedRepeatThreshold,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const Cluster &cluster,
        bool withGaps,
        FragmentMetadataList &fragments) const;

    template <typename MatchFinderT, typename FragmentCallbackT>
    AlignmentType buildAllHeadAnchored(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const unsigned filterContigId,
        const std::size_t seedRepeatThreshold,
        const unsigned headLengthMax,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const MatchFinderT &matchFinder,
        const Cluster &cluster,
        FragmentCallbackT callback) const;


    bool realignBadUngappedAlignments(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        FragmentMetadataList &fragments) const
    {
        return gappedAligner_.realignBadUngappedAlignments(
            gappedMismatchesMax_, smitWatermanGapsMax_, contigList, readMetadata, fragments, adapterClipper, cigarBuffer_);
    }

private:
    static const unsigned READS_MAX = 2;
    static const unsigned MIN_CANDIDATES = 3;
    const unsigned repeatThreshold_;
    const unsigned gappedMismatchesMax_;
    const unsigned smitWatermanGapsMax_;
    const bool noSmithWaterman_;
    const bool smartSmithWaterman_;
    const bool splitAlignments_;

    const AlignmentCfg &alignmentCfg_;

//    static std::array<std::array<std::pair<std::atomic<std::size_t>, std::atomic<std::size_t> >, 7>, 7> matchSeedCounts_;
//    static std::array<std::atomic<std::size_t>, 7> alignedSeedCounts_;
//    static std::array<std::atomic<std::size_t>, 7> seedMatchCounts_;
//    static bool countsTraced_;
    /**
     * \brief flag per seed indicating whether the seed matches are ignored due to
     *        a high repeat match
     */
    Cigar &cigarBuffer_;
    mutable Cigar oneFragmentCigarBuffer_;

    const templateBuilder::UngappedAligner ungappedAligner_;
    mutable templateBuilder::GappedAligner gappedAligner_;

    mutable ReferenceOffsetLists fwMergeBuffers_;
    mutable ReferenceOffsetLists rvMergeBuffers_;
    mutable MatchLists matchLists_;
    struct BestMatch
    {
        BestMatch (const Match &match, const unsigned mismatches):
            match_(match), mismatches_(mismatches){}
        Match match_;
        unsigned mismatches_;

        friend std::ostream & operator <<(std::ostream &os, const BestMatch &bm)
        {
            return os << "BestMatch("  << bm.match_ << "," << bm.mismatches_ << "mm)" << std::endl;
        }
    };
    mutable std::vector<BestMatch> bestMatches_;

    /**
     ** \brief add a match, either by creating a new instance of
     ** FragmentMetadata or by updating an existing one
     **
     ** Initializes the FragmentMetadata in the list for the corresponding
     ** readIndex with contigId, orientation (reverse flag) and
     ** position. The fragment is initially located at the leftmost position of
     ** the read on the forward strand of the contig. This means that the
     ** position can be negative.
     **
     ** Note: spurious FragmentMetadata creation is avoided by checking if the
     ** last FragmentMetadata created for the read has same contigId, position
     ** and orientation.
     **/
    void addMatch(
        const flowcell::ReadMetadata &readMetadata,
        const Match &match,
        const Cluster &cluster,
        FragmentMetadataList &fragments);

    AlignmentType findBestAlignments(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const Cluster &cluster,
        const bool withGaps,
        const MatchLists &matchLists,
        const unsigned uncheckedSeeds,
        FragmentMetadataList &fragments) const;

    AlignmentType findBestMatches(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const Cluster &cluster,
        const MatchLists &matchLists,
        std::vector<BestMatch> &bestMatches) const;

    FragmentMetadata makeAlignment(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const Cluster &cluster,
        FragmentSequencingAdapterClipper &adapterClipper,
        const Match &match,
        const unsigned uncheckedSeeds,
        Cigar &cigarBuffer) const;

    bool makeBestAlignments(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const Cluster &cluster,
        const bool withGaps,
        const std::vector<BestMatch> &bestMatches,
        const unsigned uncheckedSeeds,
        FragmentSequencingAdapterClipper &adapterClipper,
        FragmentMetadataList &fragments) const;

    bool updateBestMatches(
        const Match &match,
        const unsigned mismatches,
        std::vector<BestMatch> &bestMatches) const;

    bool collectBestMatches(
        const Matches& matches, const reference::ContigList& contigList, const Read& read,
        //std::size_t &counts,
        std::vector<BestMatch>& bestMatches) const;

    static const unsigned CHECK_MATCH_GROUPS_MAX = 2; // no point to go through all match groups.
//    void updateHitStats(const common::StaticVector<int,CHECK_MATCH_GROUPS_MAX>& seedCounts, std::size_t counts[CHECK_MATCH_GROUPS_MAX]) const;
};

inline double calculateLogProbability(
    const reference::Contig &contig,
    const Read &read,
    const bool reverse,
    int64_t alignmentPosition)
{
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const std::vector<char> &quality = read.getStrandQuality(reverse);

    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    templateBuilder::AlignerBase::clipReference(contig.size(), alignmentPosition, sequenceBegin, sequenceEnd);
    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);

    return FragmentMetadata::calculateLogProbability(
        std::distance(sequenceBegin, sequenceEnd),
        contig.begin() + alignmentPosition, sequenceBegin, quality.begin() + firstMappedBaseOffset);
}

/**
 * \return true if at least one fragment was built.
 */
template <typename MatchFinderT>
AlignmentType FragmentBuilder::buildBest(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const std::size_t seedRepeatThreshold,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const MatchFinderT &matchFinder,
    const Cluster &cluster,
    const bool withGaps,
    FragmentMetadataList &fragments) const
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::build: cluster " << cluster.getId() << " " << readMetadata);
    ISAAC_ASSERT_MSG(cluster.getNonEmptyReadsCount() > readMetadata.getIndex(), "cluster geometry must match");

    fragments.clear();
    ISAAC_ASSERT_MSG(!matchLists_.empty(), "empty matches lists");
    const std::size_t uncheckedSeeds = matchFinder.findReadMatches(
        contigList, cluster, readMetadata, seedRepeatThreshold, matchLists_, fwMergeBuffers_, rvMergeBuffers_);

    const AlignmentType ret = findBestAlignments(
        contigList, readMetadata, adapterClipper, cluster, withGaps, matchLists_, uncheckedSeeds, fragments);
    if (Nm == ret && uncheckedSeeds)
    {
        return Rm;
    }

    return ret;
}

/**
 * \return true if at least one fragment was built.
 */
template <typename MatchFinderT, typename FragmentCallbackT>
AlignmentType FragmentBuilder::buildAllHeadAnchored(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    const unsigned filterContigId,
    const std::size_t seedRepeatThreshold,
    const unsigned headLengthMax,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const MatchFinderT &matchFinder,
    const Cluster &cluster,
    FragmentCallbackT callback) const
{
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "FragmentBuilder::buildAllHeadAnchored: cluster " << cluster.getId());
    ISAAC_ASSERT_MSG(cluster.getNonEmptyReadsCount() > readMetadata.getIndex(), "cluster geometry must match");

    ISAAC_ASSERT_MSG(!matchLists_.empty(), "unformatted matches lists");

    const std::size_t uncheckedSeeds = matchFinder.findHeadAnchoredReadMatches(
        contigList,
        cluster, readMetadata, filterContigId, seedRepeatThreshold, headLengthMax, matchLists_, fwMergeBuffers_, rvMergeBuffers_);

    std::size_t count = 0;
    for (unsigned supportingSeeds = matchLists_.size() - 1; 0 != supportingSeeds; --supportingSeeds)
    {
        const Matches &matches = matchLists_[supportingSeeds];
        for (const Match &match: matches)
        {
            oneFragmentCigarBuffer_.clear();
            isaac::alignment::FragmentMetadata fragment = makeAlignment(contigList, readMetadata, cluster, adapterClipper, match, uncheckedSeeds, oneFragmentCigarBuffer_);
            if (fragment.isAligned())
            {
                callback(std::move(fragment));
                ++count;
            }
        }
    }

    return !count ? (uncheckedSeeds ? Rm : Nm) : Normal;
}

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_HH
