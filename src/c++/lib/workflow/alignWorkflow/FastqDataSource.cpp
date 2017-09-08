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
 ** \file FastqDataSource.cpp
 **
 ** \brief see FastqDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "workflow/alignWorkflow/FastqDataSource.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{

void FastqBaseCallsSource::generateAllFilePaths()
{
    // prepare all paths as we will not be able to allocate strings once we start loading
    for (const unsigned lane : lanes_)
    {
        lanePaths_.resize(std::max<unsigned>(lanePaths_.size(), lane + 1));
        const boost::filesystem::path read1Path = fastqFlowcellLayout_.getLaneReadAttribute<flowcell::Layout::Fastq,
            flowcell::FastqFilePathAttributeTag>(lane, fastqFlowcellLayout_.getReadMetadataList().at(0).getNumber());
        if (1 == fastqFlowcellLayout_.getReadMetadataList().size())
        {
            lanePaths_.at(lane) = std::make_pair(read1Path, boost::filesystem::path());
        }
        else
        {
            const boost::filesystem::path read2Path = fastqFlowcellLayout_.getLaneReadAttribute<flowcell::Layout::Fastq,
                flowcell::FastqFilePathAttributeTag>(lane, fastqFlowcellLayout_.getReadMetadataList().at(1).getNumber());
            lanePaths_.at(lane) = std::make_pair(read1Path, read2Path);
        }
    }
}

FastqBaseCallsSource::FastqBaseCallsSource(
    const unsigned clustersAtATimeMax,
    const unsigned coresMax,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::Layout &fastqFlowcellLayout,
    common::ThreadVector &threads) :
        tileClustersMax_(clustersAtATimeMax),
        coresMax_(coresMax),
        barcodeMetadataList_(barcodeMetadataList),
        fastqFlowcellLayout_(fastqFlowcellLayout),
        clusterLength_(flowcell::getTotalReadLength(fastqFlowcellLayout_.getReadMetadataList()) + fastqFlowcellLayout_.getBarcodeLength() + fastqFlowcellLayout_.getReadNameLength()),
        loadedClusters_(clusterLength_),
        loadingClusters_(clusterLength_),
        lanes_(fastqFlowcellLayout.getLaneIds()),
        loadingLaneIterator_(lanes_.begin()),
        fastqLoader_(fastqFlowcellLayout_.getAttribute<flowcell::Layout::Fastq, flowcell::FastqVariableLengthOk>(), 0, threads, coresMax_)
{
    loadedClusters_.reset(clusterLength_, tileClustersMax_);
    // reserve space.
    loadingClusters_.reset(clusterLength_, tileClustersMax_);
    // indicate we have nothing loaded
    loadedClusters_.reset(clusterLength_, 0);
    loadingClusters_.reset(clusterLength_, 0);
    // prepare all paths as we will not be able to allocate strings once we start loading
    generateAllFilePaths();
    tileLoadThread_ = std::thread([this](){loadTilesThread();});
}

FastqBaseCallsSource::~FastqBaseCallsSource()
{
    terminateRequested_ = true;
    stateChangeEvent_.notify_all();
    tileLoadThread_.join();
}

void FastqBaseCallsSource::reopenFastq(
    const flowcell::Layout &fastqFlowcellLayout,
    const unsigned lane,
    io::FastqLoader &fastqLoader) const
{
    if (1 == fastqFlowcellLayout.getReadMetadataList().size())
    {
        // this will keep the current files open if the paths don't change
        fastqLoader.open(lanePaths_.at(lane).first, fastqFlowcellLayout.getAttribute<flowcell::Layout::Fastq, flowcell::FastqBaseQ0>());
    }
    else // assume paired data
    {
        // this will keep the current files open if the paths don't change
        fastqLoader.open(lanePaths_.at(lane).first, lanePaths_.at(lane).second, fastqFlowcellLayout.getAttribute<flowcell::Layout::Fastq, flowcell::FastqBaseQ0>());
    }
}

void FastqBaseCallsSource::loadNextTile()
{
    // As we don't know whether the current lane has been completely loaded or we're in the
    // middle of discovering it's tiles, just attempt to load more data for it and stop only
    // when some data is loaded for this or subsequent lane or we run out of lanes.
    unsigned clustersLoaded = 0;
    while (!clustersLoaded && lanes_.end() != loadingLaneIterator_)
    {
        reopenFastq(fastqFlowcellLayout_, *loadingLaneIterator_, fastqLoader_);
        ISAAC_THREAD_CERR<< "Resetting Fastq data for " << tileClustersMax_ << " clusters" << std::endl;
        loadingClusters_.reset(clusterLength_, tileClustersMax_);
        ISAAC_THREAD_CERR<< "Resetting Fastq data done for " << loadingClusters_.getClusterCount() << " clusters" << std::endl;
        // load clusters, return tile breakdown based on tileClustersMax_
        clustersLoaded = fastqLoader_.loadClusters(tileClustersMax_, fastqFlowcellLayout_.getReadNameLength(),
                                                   fastqFlowcellLayout_.getReadMetadataList(), loadingClusters_.cluster(0));
        ISAAC_THREAD_CERR<< "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
        if (!clustersLoaded)
        {
            // there was nothing left to load for this lane, proceed with the next one.
            ++loadingLaneIterator_;
            loadingTile_ = 0;
        }
        else
        {
            loadingTile_++;
        }
    }
    loadingClusters_.reset(clusterLength_, clustersLoaded);
}

void FastqBaseCallsSource::loadTilesThread()
{
    ISAAC_THREAD_CERR<< "Fastq load thread started" << std::endl;
    std::unique_lock<std::mutex> lock(stateMutex_);
    while (!terminateRequested_ && !noMoreData_)
    {
        if (!loadingClusters_.getClusterCount())
        {
            // we have room to load more data
            try
            {
                // make sure client thread can pick up previous tile while we're loading this one
                common::unlock_guard<std::unique_lock<std::mutex> > unlock(lock);
                loadNextTile();
            }
            catch (...)
            {
                forceTermination_ = true;
                throw;
            }
        }

        if (loadedTile_)
        {
            // last loaded tile has not been picked up
            stateChangeEvent_.wait(lock);
        }

        ISAAC_ASSERT_MSG(!loadedTile_, "Expected client to be ready for next tile");

        if (lanes_.end() != loadingLaneIterator_)
        {
            ISAAC_ASSERT_MSG(!loadedClusters_.getClusterCount(), "Expected client buffer to be empty");
            loadedClusters_.swap(loadingClusters_);
            loadedTile_ = loadingTile_;
            loadedLane_ = *loadingLaneIterator_;
        }
        else
        {
            noMoreData_ = true;
        }
        stateChangeEvent_.notify_all();
    }
    ISAAC_THREAD_CERR<< "Fastq load thread terminated" << std::endl;
}

flowcell::TileMetadataList FastqBaseCallsSource::discoverTiles()
{
    std::unique_lock<std::mutex> lock(stateMutex_);

    while (!loadedTile_ && !noMoreData_)
    {
        if (forceTermination_)
        {
            BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
        }
        stateChangeEvent_.wait(lock);
    }
    flowcell::TileMetadataList ret;
    if (loadedTile_ && loadedClusters_.getClusterCount())
    {
        ret.push_back(flowcell::TileMetadata(
            fastqFlowcellLayout_.getFlowcellId(), fastqFlowcellLayout_.getIndex(), loadedTile_,
            loadedLane_, loadedClusters_.getClusterCount(), 0));
    }
    return ret;
}

void FastqBaseCallsSource::loadClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    std::unique_lock<std::mutex> lock(stateMutex_);
    if (forceTermination_)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }

    ISAAC_ASSERT_MSG(loadedTile_, "Expected to have a tile loaded by now");
    ISAAC_ASSERT_MSG(tileMetadata.getFlowcellIndex() == fastqFlowcellLayout_.getIndex(), "Unexpected tile requested");
    ISAAC_ASSERT_MSG(tileMetadata.getLane() == loadedLane_, "Unexpected tile lane requested: " << tileMetadata);
    ISAAC_ASSERT_MSG(tileMetadata.getTile() == loadedTile_, "Unexpected tile tile requested: " << tileMetadata);

    bclData.swap(loadedClusters_);
    loadedClusters_.reset(clusterLength_, 0);
    loadedTile_ = 0;
    stateChangeEvent_.notify_all();
}

} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
