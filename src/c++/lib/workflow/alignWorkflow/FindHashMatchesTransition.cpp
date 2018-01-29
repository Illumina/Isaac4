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
 ** \file FindHashMatchesTransition.cpp
 **
 ** \brief see FindHashMatchesTransition.hh
 **
 ** \author Roman Petrovski
 **/

#include <boost/ref.hpp>

#include "alignment/HashMatchFinder.hh"
#include "alignment/matchSelector/BinningFragmentStorage.hh"
#include "alignment/matchSelector/DebugStorage.hh"
#include "build/Build.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/Numa.hh"
#include "demultiplexing/DemultiplexingStatsXml.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/BclBgzfDataSource.hh"
#include "workflow/alignWorkflow/BclDataSource.hh"
#include "workflow/alignWorkflow/MultiTileDataSource.hh"
#include "workflow/alignWorkflow/FastqDataSource.hh"
#include "workflow/alignWorkflow/FindHashMatchesTransition.hh"


namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{
namespace findHashMatchesTransition
{

void wait(
    bool &signal,
    boost::condition_variable &condition,
    boost::unique_lock<boost::mutex> &lock,
    bool &forceTermination)
{
    if (forceTermination)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }

    while (signal)
    {
        if (forceTermination)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }

        condition.wait(lock);
    }
    signal = true;
}

void release(
    bool &signal,
    boost::condition_variable &condition,
    bool &forceTermination,
    const bool exceptionUnwinding)
{
    if (exceptionUnwinding)
    {
        forceTermination = true;
    }
    signal = false;
    condition.notify_all();
}

void IoOverlapThreadWorker::binQscores(alignment::BclClusters &bclData) const
{
    ISAAC_THREAD_CERR << "Binning qscores" << std::endl;

    for(std::size_t clusterId = 0; bclData.getClusterCount() != clusterId; ++clusterId)
    {
        alignment::BclClusters::iterator clusterBegin = bclData.cluster(clusterId);
        for (const flowcell::ReadMetadata &readMetadata : flowcellLayout_.getReadMetadataList())
        {
            BclClusterFields::IteratorPair pair = bclFields_.getBcl(clusterBegin, readMetadata.getIndex());
            std::for_each(
                pair.first, pair.second, [this](char &bcl){bcl = fullBclQScoreTable_[static_cast<unsigned char>(bcl)];});
        }
    }
    ISAAC_THREAD_CERR << "Binning qscores done" << std::endl;
}

/**
 * \brief Finds matches for the lane. Updates foundMatches with match information and tile metadata identified during
 *        the processing.
 */
template <typename HashMatchFinder, typename DataSourceT>
void IoOverlapThreadWorker::run(
    const flowcell::TileMetadataList &unprocessedTiles,
    unsigned &nextTile,
    unsigned &nextUnprocessedTile,
    DataSourceT &dataSource,
    const HashMatchFinder &matchFinder,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    common::ScopedMallocBlock &mallocBlock)
{
    boost::unique_lock<boost::mutex> lock(mutex_);
    while(unprocessedTiles.size() != nextTile)
    {
        const unsigned ourTile = nextTile++;
        const flowcell::TileMetadata &tileMetadata = unprocessedTiles.at(ourTile);

        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&release, boost::ref(loading_), boost::ref(stateChangedCondition_), boost::ref(forceTermination_), _1))
//        ISAAC_BLOCK_WITH_CLENAUP([this](bool){release(loading_, stateChangedCondition_);})
        {
            wait(loading_, stateChangedCondition_, lock, forceTermination_);
            {
                common::ScopedMallocBlockUnblock unblockMalloc(mallocBlock);
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);

                dataSource.resetBclData(tileMetadata, tileClusters_);
                dataSource.loadClusters(tileMetadata, tileClusters_);
                if(qScoreBin_)
                {
                    binQscores(tileClusters_);
                }
            }
        }

        ISAAC_BLOCK_WITH_CLENAUP([&](bool exceptionUnwinding)
        {
            if (exceptionUnwinding) {forceTermination_ = true;}
            stateChangedCondition_.notify_all();
        })
        {
            // make sure the order in which tiles are processed is same between different runs.
            while (nextUnprocessedTile != ourTile)
            {
                if (forceTermination_)
                {
                    BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
                }

                stateChangedCondition_.wait(lock);
            }
#ifdef ISAAC_ALIGNMENT_LOOP_ENABLED
            while (true)
#endif //ISAAC_ALIGNMENT_LOOP_ENABLED
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                matchSelector_.parallelSelect(tileClusterInfo, barcodeTemplateLengthStatistics, tileMetadata, matchFinder, tileClusters_, fragmentStorage_);
            }

            // swap the flush buffers while we still have compute lock
            wait(flushing_, stateChangedCondition_, lock, forceTermination_);
            {
                fragmentStorage_.prepareFlush();
            }
            ++nextUnprocessedTile;
        }

        // flush asynchronously so that other guys can load and align at the same time
        ISAAC_BLOCK_WITH_CLENAUP(boost::bind(&release, boost::ref(flushing_), boost::ref(stateChangedCondition_), boost::ref(forceTermination_), _1))
//        ISAAC_BLOCK_WITH_CLENAUP([&](bool){release(flushing_, stateChangedCondition_);})
        {
            // flush slot already acquired when we had the compute slot but the state could have changed in between.
            if (forceTermination_)
            {
                BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
            }
            {
                common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
                fragmentStorage_.flush();
            }
        }
    }
}

} // namespace findHashMatchesTransition

FindHashMatchesTransition::FindHashMatchesTransition(
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
    const unsigned int maxThreadCount,
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
    const alignment::AlignmentCfg &alignmentCfg,
    const alignment::TemplateBuilder::DodgyAlignmentScore dodgyAlignmentScore,
    const unsigned anomalousPairHandicap,
    const bool qScoreBin,
    const boost::array<char, 256> &fullBclQScoreTable,
    const unsigned targetBinLength,
    const uint64_t targetBinSize,
    const double expectedBgzfCompressionRatio,
    const bool preSortBins,
    const bool preAllocateBins,
    const std::string &binRegexString,
    const unsigned detectTemplateBlockSize
    )
    : hashTableBucketCount_(hashTableBucketCount)
    , flowcellLayoutList_(flowcellLayoutList)
    , tempDirectory_(tempDirectory)
    , demultiplexingStatsXmlPath_(demultiplexingStatsXmlPath)
    , coresMax_(maxThreadCount)
    , candidateMatchesMax_(candidateMatchesMax)
    , matchFinderMaxRepeats_(std::max(matchFinderTooManyRepeats, std::max(matchFinderWayTooManyRepeats, matchFinderShadowSplitRepeats)))
    , seedBaseQualityMin_(seedBaseQualityMin)
    , seedLength_(seedLength)
    , repeatThreshold_(repeatThreshold)
    , neighborhoodSizeThreshold_(neighborhoodSizeThreshold)
    , barcodeMetadataList_(barcodeMetadataList)
    , cleanupIntermediary_(cleanupIntermediary)
    , bclTilesPerChunk_(bclTilesPerChunk)
    , ignoreMissingBcls_(ignoreMissingBcls)
    , ignoreMissingFilters_(ignoreMissingFilters)
    , availableMemory_(availableMemory)
    , clustersAtATimeMax_(clustersAtATimeMax)
    , ignoreNeighbors_(ignoreNeighbors)
    , ignoreRepeats_(ignoreRepeats)
    , inputLoadersMax_(inputLoadersMax)
    , tempSaversMax_(tempSaversMax)
    , memoryControl_(memoryControl)
    , clusterIdList_(clusterIdList)
    , sortedReferenceMetadataList_(sortedReferenceMetadataList)
    , extractClusterXy_(extractClusterXy)
    , debugAlignments_(true)
    , maxGapsPerRead_(std::max<unsigned>(smitWatermanGapsMax, splitAlignments))
    , targetBinLength_(targetBinLength)
    , targetBinSize_(targetBinSize)
    , expectedBgzfCompressionRatio_(expectedBgzfCompressionRatio)
    , keepUnaligned_(keepUnaligned)
    , preSortBins_(preSortBins)
    , preAllocateBins_(preAllocateBins)
    , binRegexString_(binRegexString)

    // Have thread pool for the maximum number of threads we may potentially need.
    , threads_(std::max(inputLoadersMax_, coresMax_))

    , ioOverlapThreads_(2)
    , contigLists_(contigLists)

    , alignmentCfg_(alignmentCfg)
    , matchSelector_(
        coresMax_,
        barcodeMetadataList_,
        flowcellLayoutList_,
        contigLists_,
        candidateMatchesMax_,
        repeatThreshold_,
        clusterIdList_,
        mateDriftRange,
        userTemplateLengthStatistics,
        mapqThreshold,
        perTileTls,
        pfOnly,
        seedLength_,
        SEEDS_PER_MATCH_MAX,
        matchFinderTooManyRepeats,
        matchFinderWayTooManyRepeats,
        matchFinderShadowSplitRepeats,
        collectCycleStats,
        baseQualityCutoff,
        keepUnaligned,
        clipSemialigned,
        clipOverlapping,
        scatterRepeats,
        rescueShadows,
        trimPEAdapters,
        anchorMate,
        gappedMismatchesMax,
        smitWatermanGapsMax,
        smartSmithWaterman,
        smithWatermanGapSizeMax,
        splitAlignments,
        alignmentCfg_,
        dodgyAlignmentScore,
        anomalousPairHandicap,
        common::ScopedMallocBlock::Strict == memoryControl_,
        detectTemplateBlockSize),
        qScoreBin_(qScoreBin),
        fullBclQScoreTable_(fullBclQScoreTable)
{
}


template <typename KmerT>
void FindHashMatchesTransition::perform(
    alignWorkflow::FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath)
{
    align<KmerT>(foundMatches, binMetadataList, barcodeTemplateLengthStatistics, matchSelectorStatsXmlPath);
}


template<class It,class End>
void FindHashMatchesTransition::perform(
    const unsigned seedLength,
    alignWorkflow::FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath,
    boost::mpl::true_ endofvec)
{
    ISAAC_ASSERT_MSG(false, "Unexpected seed length " << seedLength);
}

template<class It,class End>
void FindHashMatchesTransition::perform(
    const unsigned seedLength,
    alignWorkflow::FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath,
    boost::mpl::false_)
{
    if(seedLength == boost::mpl::deref<It>::type::value)
    {
        perform<oligo::BasicKmerType<boost::mpl::deref<It>::type::value> >(
            foundMatches, binMetadataList, barcodeTemplateLengthStatistics, matchSelectorStatsXmlPath);
    }
    else
    {
        typedef typename boost::mpl::next<It>::type Next;
        perform<Next,End>(
            seedLength, foundMatches, binMetadataList, barcodeTemplateLengthStatistics, matchSelectorStatsXmlPath, typename boost::is_same<Next,End>::type());
    }
}

void FindHashMatchesTransition::perform(
    const unsigned seedLength,
    alignWorkflow::FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath)
{
    typedef boost::mpl::begin<oligo::SUPPORTED_KMERS>::type begin;
    typedef boost::mpl::end<oligo::SUPPORTED_KMERS>::type end;

    perform<begin,end>(
        seedLength, foundMatches, binMetadataList, barcodeTemplateLengthStatistics, matchSelectorStatsXmlPath,
        boost::is_same<begin,end>::type());
}


void FindHashMatchesTransition::resolveBarcodes(
    const flowcell::Layout &flowcell,
    const flowcell::BarcodeMetadataList &barcodeGroup,
    BarcodeSource &barcodeSource,
    // this contains tiles we are processing but they are not placed at the tile.getIndex()
    flowcell::TileMetadataList unprocessedTiles,
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    demultiplexing::DemultiplexingStats &demultiplexingStats)
{
    ISAAC_ASSERT_MSG(!barcodeGroup.empty(), "At least 'none' barcode must be defined");
    if (1 == barcodeGroup.size())
    {
        const flowcell::BarcodeMetadata barcode = barcodeGroup.at(0);
        ISAAC_ASSERT_MSG(barcode.isNoIndex(), "If barcode group has only one entry it must be the 'NoIndex' barcode");
        BOOST_FOREACH(const flowcell::TileMetadata &tile, unprocessedTiles)
        {
            for (unsigned clusterId = 0; clusterId < tile.getClusterCount(); ++clusterId)
            {
                tileClusterInfo.setBarcodeIndex(tile.getIndex(), clusterId, barcode.getIndex());
                demultiplexingStats.recordBarcode(demultiplexing::BarcodeId(tile.getIndex(), barcode.getIndex(), clusterId, 0));
            }
            ISAAC_THREAD_CERR << "Forced barcode index for clusters of " << tile << " to " << barcode << std::endl;
        }
    }
    else
    {
        demultiplexing::BarcodeResolver barcodeResolver(barcodeMetadataList_, barcodeGroup);

        flowcell::TileMetadataList currentTiles; currentTiles.reserve(unprocessedTiles.size());

        while (!unprocessedTiles.empty())
        {
            currentTiles.clear();
            if (!demultiplexing::BarcodeMemoryManager::selectTiles(unprocessedTiles, currentTiles))
            {
                BOOST_THROW_EXCEPTION(common::MemoryException("Insufficient memory to load barcodes even for just one tile: " +
                    boost::lexical_cast<std::string>(unprocessedTiles.back())));
            }

            demultiplexing::Barcodes barcodes;
            // this will take at most the same amount of ram as a set of singleseeds
            ISAAC_ASSERT_MSG(barcodeGroup.size(), "Barcode list must be not empty");
            ISAAC_ASSERT_MSG(barcodeGroup.at(0).isDefault(), "The very first barcode must be the 'unknown indexes or no index' one");
            barcodeSource.loadBarcodes(flowcell, barcodeGroup.at(0).getIndex(), currentTiles, barcodes);
            barcodeResolver.resolve(barcodes, demultiplexingStats);

            BOOST_FOREACH(const demultiplexing::Barcode &barcode, barcodes)
            {
                tileClusterInfo.setBarcodeIndex(barcode.getTile(), barcode.getCluster(), barcode.getBarcode());
            }
        }
    }
}

inline bool isPrime(const std::size_t num)
{
    if (num <= 3)
    {
        return num > 1;
    }
    else if (num % 2 == 0 || num % 3 == 0)
    {
        return false;
    }

    for (std::size_t i = 5; i * i <= num; i += 6)
    {
        if (num % i == 0 || num % (i + 2) == 0)
        {
            return false;
        }
    }
    return true;
}
//
//static std::size_t nextprime(std::size_t prime)
//{
//    do
//    {
//        ++prime;
//    }
//    while (!isPrime(prime));
//    return prime;
//}
//
//static std::size_t prevprime(std::size_t prime)
//{
//    do
//    {
//        --prime;
//    }
//    while (!isPrime(prime));
//    return prime;
//}

template <typename ReferenceHashT>
ReferenceHashT buildReferenceHash(
    const reference::ContigList &contigList,
    const std::size_t hashTableBucketCount,
    common::ThreadVector &threads,
    const unsigned coresMax)
{
    reference::ReferenceHasher<ReferenceHashT> hasher(contigList, threads, coresMax);

    ReferenceHashT ret = hasher.generate(hashTableBucketCount);

    return ret;
}

/**
 * \brief Finds matches for the lane. Updates foundMatches with match information and tile metadata identified during
 *        the processing.
 */
template <typename ReferenceHashT, typename DataSourceT>
void FindHashMatchesTransition::findLaneMatches(
    const ReferenceHashT &referenceHash,
    const flowcell::Layout &flowcell,
    const unsigned lane,
    const flowcell::BarcodeMetadataList &laneBarcodes,
    const flowcell::TileMetadataList &unprocessedTiles,
    DataSourceT &dataSource,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    std::vector<findHashMatchesTransition::IoOverlapThreadWorker> &ioOverlapThreadWorkers)
{
    if (!unprocessedTiles.empty())
    {
        alignment::matchFinder::TileClusterInfo tileClusterInfo(unprocessedTiles);
        ISAAC_THREAD_CERR << "Resolving barcodes for " << flowcell << " lane " << lane << std::endl;
        resolveBarcodes(
            flowcell, laneBarcodes,
            dataSource,
            unprocessedTiles, tileClusterInfo, demultiplexingStats);
        ISAAC_THREAD_CERR << "Resolving barcodes done for " << flowcell << " lane " << lane << std::endl;

        ISAAC_THREAD_CERR << "Finding hash matches with repeat threshold: " << repeatThreshold_ << std::endl;

        alignment::ClusterHashMatchFinder<ReferenceHashT, SEEDS_PER_MATCH_MAX> matchFinder(
            referenceHash, candidateMatchesMax_, seedBaseQualityMin_, matchFinderMaxRepeats_);

        matchSelector_.reserveMemory(unprocessedTiles);

        {
            unsigned current = 0;
            unsigned nextUnprocessed = 0;
            common::ScopedMallocBlock  mallocBlock(memoryControl_);
            ioOverlapThreads_.execute
            (
                [&](const unsigned threadNumber, const unsigned threadsTotal)
                {
                    ioOverlapThreadWorkers.at(threadNumber).run(
                        unprocessedTiles, current, nextUnprocessed, dataSource, matchFinder, tileClusterInfo,
                        barcodeTemplateLengthStatistics, mallocBlock);
                },
                std::min(unprocessedTiles.size(), ioOverlapThreadWorkers.size())
            );
        }

        ISAAC_THREAD_CERR << "Finding Single-seed matches done" << std::endl;
    }
}

/**
 * \brief Collect all barcodes belonging to the flowcell lane. Default barcode is placed at the beginning of the result list
 *
 * \return Returns subset of barcodeMetadataList or an empty list if none of the barcodes are mapped to a reference
 */
static flowcell::BarcodeMetadataList findFlowcellLaneBarcodes(
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout& flowcell,
    const unsigned lane)
{
    bool allBarcodesUnmapped = true;
    // put a placeholder for the unknown barcode in the beginning of the list
    flowcell::BarcodeMetadataList ret(1);
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        if (barcode.getFlowcellId() == flowcell.getFlowcellId() &&
            barcode.getLane() == lane)
        {
            ISAAC_THREAD_CERR << "adding " << barcode << std::endl;
            if (barcode.isDefault())
            {
                ISAAC_ASSERT_MSG(flowcell::BarcodeMetadata::INVALID_INDEX == ret.at(0).getIndex(), "More than one explicit specification for 'default' barcode within the group.");
                ret[0] = barcode;
            }
            else
            {
                ret.push_back(barcode);
            }
            allBarcodesUnmapped &= barcode.isUnmappedReference();
        }
    }

    ISAAC_ASSERT_MSG(flowcell::BarcodeMetadata::INVALID_INDEX != ret.at(0).getIndex(), "Missing default barcode specification");

    if (allBarcodesUnmapped)
    {
        ret.clear();
    }
    return ret;
}


template <typename ReferenceHashT, typename DataSourceT>
void FindHashMatchesTransition::processFlowcellTiles(
    const ReferenceHashT &referenceHash,
    const flowcell::Layout& flowcell,
    DataSourceT &dataSource,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    FoundMatchesMetadata &foundMatches,
    alignment::matchSelector::FragmentStorage &fragmentStorage)
{
    ISAAC_TRACE_STAT("FindHashMatchesTransition::processFlowcellTiles before fragmentStorage.reserve")
    fragmentStorage.reserve(dataSource.getMaxTileClusters());
    ISAAC_TRACE_STAT("FindHashMatchesTransition::processFlowcellTiles after fragmentStorage.reserve")

    loading_ = false;
    flushing_ = false;
    forceTermination_ = false;

    ISAAC_TRACE_STAT("FindHashMatchesTransition::findLaneMatches before threads allocation")
    std::vector<findHashMatchesTransition::IoOverlapThreadWorker> ioOverlapThreadWorkers(
        ioOverlapThreads_.size(),
        findHashMatchesTransition::IoOverlapThreadWorker(
            mutex_,
            stateChangedCondition_,
            loading_,
            flushing_,
            forceTermination_,
            flowcell,
            dataSource.getMaxTileClusters(),
            DataSourceTraits<DataSourceT>::SUPPORTS_XY && extractClusterXy_,
            qScoreBin_,
            fullBclQScoreTable_,
            matchSelector_,
            fragmentStorage));

    ISAAC_TRACE_STAT("FindHashMatchesTransition::findLaneMatches after allocation")


    for (flowcell::TileMetadataList laneTiles = dataSource.discoverTiles(); !laneTiles.empty();
        laneTiles = dataSource.discoverTiles())
    {
        const unsigned lane = laneTiles[0].getLane();
        flowcell::BarcodeMetadataList laneBarcodes = findFlowcellLaneBarcodes(barcodeMetadataList_, flowcell, lane);
        if (laneBarcodes.empty())
        {
            ISAAC_THREAD_CERR << "Skipping flowcell " << flowcell.getFlowcellId() << " lane " << lane << " as none of the barcodes map to the reference" << std::endl;
        }
        else
        {
            BOOST_FOREACH(TileMetadata &tileMetadata, laneTiles)
            {
                foundMatches.addTile(tileMetadata);
                // this fixes the tile index to be correct in the context of the global tile list.
                // it is important for findLaneMatches to use the global tile index so that it can store
                // statistics properly
                tileMetadata = foundMatches.tileMetadataList_.back();
            }
            ISAAC_TRACE_STAT("FindHashMatchesTransition::processFlowcellTiles before findLaneMatches")
            findLaneMatches(
                referenceHash, flowcell, lane, laneBarcodes, laneTiles, dataSource,
                demultiplexingStats, barcodeTemplateLengthStatistics, ioOverlapThreadWorkers);
        }
    }
}

template <typename ReferenceHashT>
void FindHashMatchesTransition::alignFlowcells(
    const ReferenceHashT &referenceHash,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &foundMatches,
    alignment::matchSelector::FragmentStorage &fragmentStorage)
{

    BOOST_FOREACH(const flowcell::Layout& flowcell, flowcellLayoutList_)
    {
        switch (flowcell.getFormat())
        {
            case flowcell::Layout::Bam:
            {
                BackgroundBamBaseCallsSource dataSource(
                    tempDirectory_,
                    availableMemory_,
                    clustersAtATimeMax_,
                    cleanupIntermediary_,
                    // the loading itself occurs on one thread at a time only. So, the real limit is to avoid using
                    // more cores for decompression than the system actually has.
                    // On the other hand, there might be a need to limit the io to 1 thread, while allowing
                    // for the multithreaded processing of other cpu-demanding things.
                    std::min(inputLoadersMax_, coresMax_),
                    flowcell, threads_);
                processFlowcellTiles(referenceHash, flowcell, dataSource, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            case flowcell::Layout::Fastq:
            {
                FastqBaseCallsSource dataSource(
                    clustersAtATimeMax_,
                    coresMax_,
                    barcodeMetadataList_,
                    flowcell,
                    threads_);

                processFlowcellTiles(referenceHash, flowcell, dataSource, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            case flowcell::Layout::Bcl:
            {
                ISAAC_TRACE_STAT("FindHashMatchesTransition::alignFlowcells before BclBaseCallsSource()")
                BclBaseCallsSource baseCalls(
                    flowcell, ignoreMissingBcls_, ignoreMissingFilters_, threads_, inputLoadersMax_, extractClusterXy_);

                MultiTileBaseCallsSource<BclBaseCallsSource> multitileBaseCalls(bclTilesPerChunk_, flowcell, baseCalls);

                processFlowcellTiles(
                    referenceHash, flowcell, multitileBaseCalls, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            case flowcell::Layout::BclBgzf:
            {
                BclBgzfBaseCallsSource baseCalls(
                    flowcell, ignoreMissingBcls_, ignoreMissingFilters_, threads_, inputLoadersMax_, extractClusterXy_);
                MultiTileBaseCallsSource<BclBgzfBaseCallsSource> multitileBaseCalls(
                    bclTilesPerChunk_, flowcell, baseCalls);

                processFlowcellTiles(referenceHash, flowcell, multitileBaseCalls, demultiplexingStats, barcodeTemplateLengthStatistics, foundMatches, fragmentStorage);
                break;
            }

            default:
            {
                ISAAC_ASSERT_MSG(false, "Unexpected flowcell format " << flowcell.getFormat());
                break;
            }
        }
    }
}
template <typename ReferenceHashT>
void FindHashMatchesTransition::alignFlowcells(
    const ReferenceHashT &referenceHash,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    demultiplexing::DemultiplexingStats &demultiplexingStats,
    FoundMatchesMetadata &ret)
{
    // unit of genome to use for counting alignment distribution
    static const unsigned TRACKING_BIN_LENGTH = 10000;
    alignment::matchSelector::BinIndexMap binIndexMap(sortedReferenceMetadataList_.front(), TRACKING_BIN_LENGTH);

    ISAAC_TRACE_STAT("AlignWorkflow::selectMatches ")
    ISAAC_THREAD_CERR << "Selecting matches using " << binIndexMap << std::endl;

    alignment::matchSelector::BinningFragmentStorage fragmentStorage(
        tempDirectory_, keepUnaligned_, binIndexMap, sortedReferenceMetadataList_.front().getContigs(),
        barcodeMetadataList_, preAllocateBins_, targetBinSize_, targetBinLength_,
        coresMax_, binMetadataList);

#ifdef ISAAC_DEV_STATS_ENABLED
        alignment::matchSelector::DebugStorage debugStorage(
            contigLists_.node0Container().front(),
            alignmentCfg_, flowcellLayoutList_, demultiplexingStatsXmlPath_.parent_path(), barcodeMetadataList_,
            coresMax_, fragmentStorage);
        alignFlowcells(referenceHash, barcodeTemplateLengthStatistics, demultiplexingStats, ret, debugStorage);
        debugStorage.close();
#else
        alignFlowcells(referenceHash, barcodeTemplateLengthStatistics, demultiplexingStats, ret, fragmentStorage);
        fragmentStorage.close();
#endif

}

template <typename KmerT>
void FindHashMatchesTransition::align(
    FoundMatchesMetadata &foundMatches,
    alignment::BinMetadataList &binMetadataList,
    std::vector<alignment::TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const boost::filesystem::path &matchSelectorStatsXmlPath)
{
//    typedef reference::ReferenceHash<KmerT, common::NumaAllocator<void, 0> > ReferenceHash;
//    typedef reference::NumaReferenceHash<ReferenceHash> NumaReferenceHash;
//    const NumaReferenceHash referenceHash(buildReferenceHash<ReferenceHash>(contigLists_.node0Container().front(), threads_, coresMax_));

    typedef reference::ReferenceHash<KmerT, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > ReferenceHash;
    const ReferenceHash referenceHash(buildReferenceHash<ReferenceHash>(
        contigLists_.node0Container().front(), hashTableBucketCount_, threads_, coresMax_));

    FoundMatchesMetadata ret(tempDirectory_, barcodeMetadataList_, 1, sortedReferenceMetadataList_);
    demultiplexing::DemultiplexingStats demultiplexingStats(flowcellLayoutList_, barcodeMetadataList_);

    alignFlowcells(
        referenceHash, binMetadataList,
        barcodeTemplateLengthStatistics, demultiplexingStats, ret);

    dumpStats(demultiplexingStats, ret.tileMetadataList_);
    foundMatches.swap(ret);

    matchSelector_.unreserve();

    matchSelector_.dumpStats(matchSelectorStatsXmlPath);
}

void FindHashMatchesTransition::dumpStats(
    const demultiplexing::DemultiplexingStats &demultiplexingStats,
    const flowcell::TileMetadataList &tileMetadataList) const
{
    demultiplexing::DemultiplexingStatsXml statsXml;

    const unsigned maxLaneNumber = flowcell::getMaxLaneNumber(flowcellLayoutList_);
    for (unsigned lane = 1; lane <= maxLaneNumber; ++lane)
    {
        typedef std::map<std::string, demultiplexing::LaneBarcodeStats> SampleLaneBarcodeStats;
        typedef std::map<std::string, SampleLaneBarcodeStats> ProjectSampleLaneBarcodeStats;
        typedef std::map<std::string, ProjectSampleLaneBarcodeStats> FlowcellProjectSampleLaneBarcodeStats;
        FlowcellProjectSampleLaneBarcodeStats flowcellProjectSampleStats;
        BOOST_FOREACH(const flowcell::BarcodeMetadata& barcode, barcodeMetadataList_)
        {
            if (barcode.getLane() == lane)
            {
                // put one lane stat for each unknown barcode found.
                if (barcode.isUnknown())
                {
                    const flowcell::Layout& flowcell = flowcellLayoutList_.at(barcode.getFlowcellIndex());
                    statsXml.addFlowcellLane(flowcell, lane,
                                             demultiplexingStats.getLaneUnknwonBarcodeStat(barcode.getIndex()));
                }
                const demultiplexing::LaneBarcodeStats &stat = demultiplexingStats.getLaneBarcodeStat(barcode);
                flowcellProjectSampleStats[barcode.getFlowcellId()][barcode.getProject()][barcode.getSampleName()] += stat;
                flowcellProjectSampleStats[barcode.getFlowcellId()][barcode.getProject()]["all"] += stat;
                flowcellProjectSampleStats[barcode.getFlowcellId()]["all"]["all"] += stat;
                flowcellProjectSampleStats["all"][barcode.getProject()][barcode.getSampleName()] += stat;
                flowcellProjectSampleStats["all"][barcode.getProject()]["all"] += stat;
                flowcellProjectSampleStats["all"]["all"]["all"] += stat;
                statsXml.addLaneBarcode(barcode.getFlowcellId(), barcode.getProject(), barcode.getSampleName(), barcode.getName(), lane, stat);
            }
        }
        BOOST_FOREACH(const FlowcellProjectSampleLaneBarcodeStats::value_type &flowcellStats, flowcellProjectSampleStats)
        {
            BOOST_FOREACH(const ProjectSampleLaneBarcodeStats::value_type &projectStats, flowcellStats.second)
            {
                BOOST_FOREACH(const SampleLaneBarcodeStats::value_type &sampleStats, projectStats.second)
                {
                    statsXml.addLaneBarcode(flowcellStats.first, projectStats.first, sampleStats.first, "all", lane, sampleStats.second);
                }
            }
        }
    }

    std::ofstream os(demultiplexingStatsXmlPath_.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + demultiplexingStatsXmlPath_.string()));
    }
    if (!(os << statsXml)) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: failed to store MatchFinder statistics in : " + demultiplexingStatsXmlPath_.string()));
    }
}

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
