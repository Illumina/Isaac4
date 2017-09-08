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
 ** \file MatchSelector.hh
 **
 ** \brief Selection the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_HH

#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "alignment/BclClusters.hh"
#include "alignment/Cluster.hh"
#include "alignment/templateBuilder/FragmentBuilder.hh"
#include "alignment/TemplateBuilder.hh"
#include "alignment/Match.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/matchSelector/MatchSelectorStats.hh"
#include "alignment/matchSelector/SemialignedEndsClipper.hh"
#include "alignment/matchSelector/OverlappingEndsClipper.hh"
#include "alignment/matchSelector/TemplateDetector.hh"
#include "common/Threads.hpp"
#include "flowcell/BarcodeMetadata.hh"
#include "reference/Contig.hh"

namespace isaac
{
namespace alignment
{

namespace bfs = boost::filesystem;

class MatchSelector: boost::noncopyable
{
public:
    typedef flowcell::TileMetadataList TileMetadataList;
    typedef flowcell::ReadMetadataList ReadMetadataList;
    /// Construction of an instance for a given reference
    MatchSelector(
        unsigned int maxThreadCount,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const reference::NumaContigLists &contigLists,
        const std::size_t candidateMatchesMax,
        const unsigned repeatThreshold,
        const std::vector<std::size_t> &clusterIdList,
        const int mateDriftRange,
        const TemplateLengthStatistics &defaultTemplateLengthStatistics,
        const int mapqThreshold,
        const bool perTileTls,
        const bool pfOnly,
        const unsigned seedLength,
        const unsigned maxSeedsPerMatch,
        const unsigned matchFinderTooManyRepeats,
        const unsigned matchFinderWayTooManyRepeats,
        const unsigned matchFinderShadowSplitRepeats,
        const bool collectCycleStats,
        const unsigned baseQualityCutoff,
        const bool keepUnaligned,
        const bool clipSemialigned,
        const bool clipOverlapping,
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
        const TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
        const unsigned anomalousPairHandicap,
        const bool reserveBuffers,
        const unsigned detectTemplateBlockSize);

    /**
     * \brief frees the major memory reservations to make it safe to use dynamic memory allocations again
     */
    void unreserve()
    {
        threadTemplateBuilders_.clear();
        std::vector<Cluster>().swap(threadCluster_);
        std::vector<matchSelector::MatchSelectorStats>().swap(threadStats_);
    }

    void dumpStats(const boost::filesystem::path &statsXmlPath);
    void reserveMemory(
        const flowcell::TileMetadataList &tileMetadataList);

    template <typename MatchFinderT>
    void parallelSelect(
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const flowcell::TileMetadata &tileMetadata,
        const MatchFinderT &matchFinder,
        const BclClusters &bclData,
        matchSelector::FragmentStorage &fragmentStorage);

private:
    // The threading code in selectTileMatches can not deal with exception cleanup. Let it just crash for now.
    common::UnsafeThreadVector computeThreads_;

    TileMetadataList tileMetadataList_;
    /**
     * \brief threadBclFilePaths_ gets resized for every tile total read length. If the tile read lengths
     *        changes from lower to bigger, more threadBclFilePaths_ strings get allocated which breaks the whole
     *        concept of allocating things once. For now this list contains tiles in the processing order so
     *        that the total read length goes only down. TODO: cleanup this mess for example by creating
     *        MatchSelector only for the group of tiles that have the same geometry.
     */
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const reference::NumaContigLists &contigLists_;
    const unsigned repeatThreshold_;
    const std::vector<size_t> &clusterIdList_;

    const int mapqThreshold_;
    const bool pfOnly_;
    const bool collectCycleStats_;
    const unsigned baseQualityCutoff_;
    const bool keepUnaligned_;
    const bool clipSemialigned_;
    const bool clipOverlapping_;
    const std::vector<SequencingAdapterList> barcodeSequencingAdapters_;

    std::vector<matchSelector::MatchSelectorStats> allStats_;
    std::vector<matchSelector::MatchSelectorStats> threadStats_;

    std::vector<Cluster> threadCluster_;
    boost::ptr_vector<TemplateBuilder> threadTemplateBuilders_;
    std::vector<matchSelector::SemialignedEndsClipper> threadSemialignedEndsClippers_;
    std::vector<matchSelector::OverlappingEndsClipper> threadOverlappingEndsClippers_;
    // updated for barcodes relevant for the current tile
    std::vector<RestOfGenomeCorrection> restOfGenomeCorrections_;

    matchSelector::TemplateDetector templateDetector_;

    mutable boost::mutex mutex_;

    template <typename MatchFinderT>
    void alignThread(
        const unsigned threadNumber,
        const flowcell::TileMetadata & tileMetadata,
        const matchFinder::ClusterInfos &clusterInfos,
        unsigned &clusterId,
        const MatchFinderT &matchFinder,
        const BclClusters &bclData,
        const std::vector<TemplateLengthStatistics> & templateLengthStatistics,
        matchSelector::FragmentStorage &fragmentStorage);


    /**
     * \brief Construct the contig list from the SortedReference XML
     */
    reference::ContigList getContigList(
        const reference::SortedReferenceMetadata &sortedReferenceMetadata) const;

    /**
     ** \brief Helper method to generate the 'rest of the genome' correction for
     ** uniquely aligned reads and fragments.
     **
     ** There is one value for each individual reads in the readMetadataList (at
     ** the corresponding location) and one additional value for cases where all
     ** the reads match uniquely.
     **/
    std::vector<double> getRestOfGenomeCorrectionList(
        const std::vector<flowcell::ReadMetadata> &readMetadataList) const;

    template <typename MatchFinderT>
    matchSelector::TemplateAlignmentType alignCluster(
        const reference::ContigList& barcodeContigList,
        const flowcell::ReadMetadataList& tileReads,
        const SequencingAdapterList& sequencingAdapters,
        const TemplateLengthStatistics& templateLengthStatistics,
        const uint64_t barcodeIndex,
        const MatchFinderT &matchFinder,
        const RestOfGenomeCorrection& restOfGenomeCorrection,
        const unsigned threadNumber, TemplateBuilder& templateBuilder,
        const Cluster& cluster, BamTemplate& bamTemplate,
        matchSelector::MatchSelectorStats& stats,
        matchSelector::FragmentStorage &fragmentStorage);

    static const unsigned CLUSTERS_AT_A_TIME = 10000;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_HH
