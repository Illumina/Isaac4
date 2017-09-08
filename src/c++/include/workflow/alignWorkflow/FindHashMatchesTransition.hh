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
 ** \file FindHashMatchesTransition.hh
 **
 ** \brief Top level component to control the analysis process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FIND_HASH_MATCHES_TRANSITION_HH
#define iSAAC_WORKFLOW_ALIGN_WORKFLOW_FIND_HASH_MATCHES_TRANSITION_HH

#include "alignment/MatchDistribution.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "alignment/HashMatchFinder.hh"
#include "alignment/MatchSelector.hh"
#include "common/Threads.hpp"
#include "demultiplexing/BarcodeLoader.hh"
#include "demultiplexing/BarcodeResolver.hh"
#include "demultiplexing/DemultiplexingStats.hh"
#include "flowcell/Layout.hh"
#include "flowcell/BarcodeMetadata.hh"
#include "flowcell/ReadMetadata.hh"
#include "flowcell/TileMetadata.hh"
#include "oligo/Kmer.hh"
#include "reference/ReferenceHasher.hh"
#include "reference/ReferenceMetadata.hh"
#include "reference/SortedReferenceMetadata.hh"

#include "workflow/alignWorkflow/BclDataSource.hh"
#include "workflow/alignWorkflow/DataSource.hh"
#include "workflow/alignWorkflow/FoundMatchesMetadata.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace bfs = boost::filesystem;

namespace findHashMatchesTransition
{




class IoOverlapThreadWorker
{
public:
    IoOverlapThreadWorker(
        boost::mutex &mutex,
        boost::condition_variable &stateChangedCondition,
        bool &loading,
        bool &flushing,
        bool &forceTermination,
        const flowcell::Layout &flowcellLayout,
        const unsigned maxTileClusters,
        const bool extractClusterXy,
        const bool qScoreBin,
        const boost::array<char, 256> &fullBclQScoreTable,
        alignment::MatchSelector &matchSelector,
        alignment::matchSelector::FragmentStorage &fragmentStorage):
            mutex_(mutex),
            stateChangedCondition_(stateChangedCondition),
            loading_(loading),
            flushing_(flushing),
            forceTermination_(forceTermination),
            flowcellLayout_(flowcellLayout),
            qScoreBin_(qScoreBin),
            fullBclQScoreTable_(fullBclQScoreTable),
            matchSelector_(matchSelector),
            fragmentStorage_(fragmentStorage),
            tileClusters_(flowcell::getTotalReadLength(flowcellLayout.getReadMetadataList()) + flowcellLayout.getBarcodeLength() + flowcellLayout.getReadNameLength()),
            bclFields_(flowcellLayout.getReadMetadataList(), flowcellLayout.getBarcodeLength())

    {
        tileClusters_.reserveClusters(maxTileClusters, extractClusterXy);
    }

    template <typename HashMatchFinder, typename DataSourceT>
    void run(
        const flowcell::TileMetadataList &unprocessedTiles,
        unsigned &currentTile,
        unsigned &nextUnprocessedTile,
        DataSourceT &dataSource,
        const HashMatchFinder &matchFinder,
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        common::ScopedMallocBlock &mallocBlock);
private:
    boost::mutex &mutex_;
    boost::condition_variable &stateChangedCondition_;
    bool &loading_;
    bool &flushing_;
    bool &forceTermination_;
    const flowcell::Layout &flowcellLayout_;
    const bool qScoreBin_;
    const boost::array<char, 256> &fullBclQScoreTable_;

    alignment::MatchSelector &matchSelector_;
    alignment::matchSelector::FragmentStorage &fragmentStorage_;

    alignment::BclClusters tileClusters_;
    typedef alignment::BclClusterFields<alignment::BclClusters::iterator> BclClusterFields;
    BclClusterFields bclFields_;

    void binQscores(alignment::BclClusters &bclData) const;
};
} // namespace findHashMatchesTransition

class FindHashMatchesTransition: boost::noncopyable
{
public:
    typedef flowcell::TileMetadata TileMetadata;
    typedef flowcell::ReadMetadata ReadMetadata;
    typedef flowcell::ReadMetadataList ReadMetadataList;

    FindHashMatchesTransition(
        const std::size_t hashTableBucketCount,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const bool cleanupIntermediary,
        const unsigned bclTilesPerChunk,
        const bool ignoreMissingBcls,
        const bool ignoreMissingFilters,
        const uint64_t availableMemory,
        const unsigned clustersAtATimeMax,
        const bfs::path &tempDirectory,
        const bfs::path &demultiplexingStatsXmlPath,
        const unsigned maxThreadCount,
        const unsigned seedLength,
        const std::size_t candidateMatchesMax,
        const unsigned matchFinderTooManyRepeats,
        const unsigned matchFinderWayTooManyRepeats,
        const unsigned matchFinderShadowSplitRepeats,
        const unsigned seedBaseQualityMin,
        const unsigned repeatThreshold,
        const unsigned neighborhoodSizeThreshold,
        const bool ignoreNeighbors,
        const bool ignoreRepeats,
        const unsigned inputLoadersMax,
        const unsigned tempSaversMax,
        const common::ScopedMallocBlock::Mode memoryControl,
        const std::vector<std::size_t> &clusterIdList,
        const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
        const reference::NumaContigLists &contigLists,
        const bool extractClusterXy,
        const int mateDriftRange,
        const alignment::TemplateLengthStatistics &userTemplateLengthStatistics,
        const int mapqThreshold,
        const bool perTileTls,
        const bool pfOnly,
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
        const int gapMatchScore,
        const int gapMismatchScore,
        const int gapOpenScore,
        const int gapExtendScore,
        const int minGapExtendScore,
        const unsigned splitGapLength,
        const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
        const unsigned anomalousPairScoreMin,
        const bool qScoreBin,
        const boost::array<char, 256> &fullBclQScoreTable,
        const unsigned targetBinLength,
        const uint64_t targetBinSize,
        const double expectedBgzfCompressionRatio,
        const bool preSortBins,
        const bool preAllocateBins,
        const std::string &binRegexString,
        const unsigned detectTemplateBlockSize);

    template <typename KmerT>
    void perform(
        alignWorkflow::FoundMatchesMetadata &foundMatches,
        alignment::BinMetadataList &binMetadataList,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const boost::filesystem::path &matchSelectorStatsXmlPath);

    template<class It, class End>
    void perform(
        const unsigned seedLength,
        alignWorkflow::FoundMatchesMetadata &foundMatches,
        alignment::BinMetadataList &binMetadataList,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const boost::filesystem::path &matchSelectorStatsXmlPath,
        boost::mpl::true_ endofvec);

    template<class It, class End>
    void perform(
        const unsigned seedLength,
        alignWorkflow::FoundMatchesMetadata &foundMatches,
        alignment::BinMetadataList &binMetadataList,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const boost::filesystem::path &matchSelectorStatsXmlPath,
        boost::mpl::false_);

    void perform(
        const unsigned seedLength,
        alignWorkflow::FoundMatchesMetadata &foundMatches,
        alignment::BinMetadataList &binMetadataList,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const boost::filesystem::path &matchSelectorStatsXmlPath);

private:
    template<class Archive> friend void serialize(Archive & ar, FindHashMatchesTransition &, const unsigned int file_version);

    static const unsigned SEEDS_PER_MATCH_MAX = 4;
    const std::size_t hashTableBucketCount_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const bfs::path tempDirectory_;
    const bfs::path demultiplexingStatsXmlPath_;
    const unsigned coresMax_;
    const std::size_t candidateMatchesMax_;
    const unsigned matchFinderMaxRepeats_;
    const unsigned seedBaseQualityMin_;
    const unsigned seedLength_;
    const unsigned repeatThreshold_;
    const unsigned neighborhoodSizeThreshold_;
    const flowcell::BarcodeMetadataList &barcodeMetadataList_;
    const bool cleanupIntermediary_;
    const unsigned bclTilesPerChunk_;
    const bool ignoreMissingBcls_;
    const bool ignoreMissingFilters_;
    const uint64_t availableMemory_;
    const unsigned clustersAtATimeMax_;
    const bool ignoreNeighbors_;
    const bool ignoreRepeats_;
    const unsigned inputLoadersMax_;
    const unsigned tempSaversMax_;
    const common::ScopedMallocBlock::Mode memoryControl_;
    const std::vector<size_t> &clusterIdList_;

    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList_;
    const bool extractClusterXy_;

    const bool debugAlignments_;
    const unsigned maxGapsPerRead_;
    const unsigned targetBinLength_;
    const uint64_t targetBinSize_;
    const double expectedBgzfCompressionRatio_;
    const bool keepUnaligned_;
    const bool preSortBins_;
    const bool preAllocateBins_;
    const std::string &binRegexString_;

    common::ThreadVector threads_;
    common::ThreadVector ioOverlapThreads_;

    boost::mutex mutex_;
    boost::condition_variable stateChangedCondition_;
    bool loading_ = false;
    bool flushing_ = false;
    bool forceTermination_ = false;

    const isaac::reference::NumaContigLists &contigLists_;
    alignment::AlignmentCfg alignmentCfg_;
    alignment::MatchSelector matchSelector_;
    bool qScoreBin_;
    const boost::array<char, 256> &fullBclQScoreTable_;


    template <typename KmerT>
    void align(
        FoundMatchesMetadata &foundMatches,
        alignment::BinMetadataList &binMetadataList,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        const boost::filesystem::path &matchSelectorStatsXmlPath);

    template <typename ReferenceHashT>
    void alignFlowcells(
        const ReferenceHashT &referenceHash,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        FoundMatchesMetadata &foundMatches,
        alignment::matchSelector::FragmentStorage &fragmentStorage);

    template <typename ReferenceHashT>
    void alignFlowcells(
        const ReferenceHashT &referenceHash,
        alignment::BinMetadataList &binMetadataList,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        FoundMatchesMetadata &ret);

    void resolveBarcodes(
        const flowcell::Layout &flowcell,
        const flowcell::BarcodeMetadataList &barcodeGroup,
        BarcodeSource &barcodeSource,
        flowcell::TileMetadataList unprocessedTiles,
        alignment::matchFinder::TileClusterInfo &tileClusterInfo,
        demultiplexing::DemultiplexingStats &demultiplexingStats);

    template <typename ReferenceHashT, typename DataSourceT>
    void findLaneMatches(
        const ReferenceHashT &referenceHash,
        const flowcell::Layout &flowcell,
        const unsigned lane,
        const flowcell::BarcodeMetadataList &barcodeGroup,
        const flowcell::TileMetadataList &unprocessedTiles,
        DataSourceT &dataSource,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        std::vector<findHashMatchesTransition::IoOverlapThreadWorker> &ioOverlapThreadWorkers);

    template <typename ReferenceHashT, typename DataSourceT>
    void processFlowcellTiles(
        const ReferenceHashT &referenceHash,
        const flowcell::Layout& flowcell,
        DataSourceT &dataSource,
        demultiplexing::DemultiplexingStats &demultiplexingStats,
        std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
        FoundMatchesMetadata &foundMatches,
        alignment::matchSelector::FragmentStorage &fragmentStorage);

    void dumpStats(
        const demultiplexing::DemultiplexingStats &demultiplexingStats,
        const flowcell::TileMetadataList &tileMetadataList) const;
};

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac

#endif // #ifndef iSAAC_WORKFLOW_ALIGN_WORKFLOW_FIND_HASH_MATCHES_TRANSITION_HH
