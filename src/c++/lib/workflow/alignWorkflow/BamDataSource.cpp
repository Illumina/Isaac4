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
 ** \file BamDataSource.cpp
 **
 ** \brief see BamDataSource.hh
 **
 ** \author Roman Petrovski
 **/

#include "oligo/Nucleotides.hh"
#include "workflow/alignWorkflow/BamDataSource.hh"
#include "workflow/alignWorkflow/bamDataSource/SingleEndClusterExtractor.hh"

namespace isaac
{
namespace workflow
{
namespace alignWorkflow
{


/**
 * \param clusterCount  Maximum number of clusters to load
 * \param clusterIt     Insert iterator for the buffer that is sufficient to load the clusterCount
 *                      clusters of clusterLength
 * \param pfIt          Insert iterator for the buffer that is sufficient to load the clusterCount
 *                      clusters of pf flags
 *
 * \return Actual number of loaded clusters
 */
template <typename ClusterInsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadPairedReads(
    unsigned clusterCount, const unsigned nameLengthMax,
    const flowcell::ReadMetadataList &readMetadataList,
    ClusterInsertIt &clusterIt, PfInserIt &pfIt)
{
    const unsigned requestedClusterCount = clusterCount;

    while (clusterCount)
    {
        if(clusterExtractor_.extractingUnpaired())
        {
            ISAAC_THREAD_CERR << "extracting unpaired " << std::endl;
            clusterCount = clusterExtractor_.extractUnpaired(
                readMetadataList.at(0).getLength(), readMetadataList.at(1).getLength(), nameLengthMax, clusterCount,
                clusterIt, pfIt);
            // Either no room in result buffer or no more data available.
            break;
        }
        else
        {
            if (!clusterExtractor_.isEmpty())
            {
                ISAAC_THREAD_CERR << "resuming from " << clusterExtractor_.size() << " pending elements" << std::endl;
                clusterCount = clusterExtractor_.extractPairedReads(
                    nameLengthMax, clusterCount, clusterIt, pfIt, readMetadataList);

                if (!clusterCount)
                {
                    // there is no room to extract any more pending items. Resume on next call.
                    return requestedClusterCount;
                }
            }

            bamLoader_.load
            (
                boost::make_tuple(
                    boost::bind(&bamDataSource::PairedEndClusterExtractor::append<ClusterInsertIt, PfInserIt>,
                                &clusterExtractor_, _1, _2, nameLengthMax, boost::ref(clusterCount),
                                boost::ref(readMetadataList), boost::ref(clusterIt), boost::ref(pfIt)),
                    boost::bind(&bamDataSource::PairedEndClusterExtractor::removeOld,
                                &clusterExtractor_, _1, _2, boost::ref(readMetadataList)))
            );


            if (clusterCount)
            {
                // we've ran out of data in the bam file. See if unpaired items can be paired
                clusterExtractor_.startExtractingUnpaired();
            }
        }
    }

    ISAAC_THREAD_CERR << "loadPairedReads done clusterCount:" << clusterCount << std::endl;
    return requestedClusterCount - clusterCount;
}

void BamClusterLoader::open(
    const std::string &flowcellId,
    const boost::filesystem::path &bamPath)
{
    if (flowcellId_ != flowcellId)
    {
        clusterExtractor_.open(flowcellId);
        bamLoader_.open(bamPath);
        flowcellId_ = flowcellId;
    }
    else
    {
        ISAAC_THREAD_CERR << "Keeping bam stream open for flowcellId " << flowcellId_ << " " << bamPath << std::endl;
    }
}

template <typename InsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadSingleReads(
    unsigned clusterCount, const unsigned nameLengthMax, const flowcell::ReadMetadataList &readMetadataList,
    InsertIt &clusterIt, PfInserIt &pfIt)
{
    const unsigned requestedClusterCount = clusterCount;
    bamDataSource::SingleEndClusterExtractor extractor;
    bamLoader_.load
    (
        boost::make_tuple(
            boost::bind(&bamDataSource::SingleEndClusterExtractor::extractSingleRead<InsertIt, PfInserIt>,
                &extractor, _1, nameLengthMax, boost::ref(clusterCount),
                boost::ref(readMetadataList), boost::ref(clusterIt), boost::ref(pfIt)),
            boost::bind(&bamDataSource::SingleEndClusterExtractor::nothing))
    );

    ISAAC_THREAD_CERR << "loadSingleReads done clusterCount:" << clusterCount << " bgzfReader_.isEof():" << std::endl;

    return requestedClusterCount - clusterCount;}

template <typename ClusterInsertIt, typename PfInserIt>
unsigned BamClusterLoader::loadClusters(
    unsigned clusterCount, const unsigned nameLengthMax, const flowcell::ReadMetadataList &readMetadataList,
    ClusterInsertIt &clusterIt, PfInserIt &pfIt)
{
    return 2 == readMetadataList.size() ?
        loadPairedReads(clusterCount, nameLengthMax, readMetadataList, clusterIt, pfIt) :
        loadSingleReads(clusterCount, nameLengthMax, readMetadataList, clusterIt, pfIt);
}

inline std::size_t getBamFileSize(const flowcell::Layout &flowcell)
{
    return common::getFileSize(flowcell.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>().c_str());
}

BamBaseCallsSource::BamBaseCallsSource(
    const boost::filesystem::path &tempDirectoryPath,
    const uint64_t availableMemory,
    const unsigned clustersAtATimeMax,
    const bool cleanupIntermediary,
    const unsigned coresMax,
    const flowcell::Layout &bamFlowcellLayout,
    common::ThreadVector &threads) :
        bamFlowcellLayout_(bamFlowcellLayout),
        bamPath_(bamFlowcellLayout_.getAttribute<flowcell::Layout::Bam, flowcell::BamFilePathAttributeTag>()),
        tileClustersMax_(clustersAtATimeMax),
        coresMax_(coresMax),
        clusterLength_(flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList()) + bamFlowcellLayout_.getBarcodeLength() + bamFlowcellLayout_.getReadNameLength()),
        clusters_(clusterLength_),
        currentTile_(0),
        threads_(threads),
        bamClusterLoader_(
            cleanupIntermediary, 0, threads, coresMax, tempDirectoryPath,
            getBamFileSize(bamFlowcellLayout_), bamFlowcellLayout.getFlowcellId().length(),
            bamFlowcellLayout_.getReadNameLength(),
            flowcell::getTotalReadLength(bamFlowcellLayout.getReadMetadataList()),
            flowcell::getMinReadLength(bamFlowcellLayout.getReadMetadataList()))
{
}

flowcell::TileMetadataList BamBaseCallsSource::discoverTiles()
{
    flowcell::TileMetadataList ret;
    discoverTiles(ret);
    return ret;
}

unsigned BamBaseCallsSource::loadNextTile()
{
    const unsigned clustersToLoad = tileClustersMax_;
    //    const unsigned clustersToLoad = clustersAtATimeMax_;
    clusters_.reset(clusterLength_, clustersToLoad);
    // load clusters, return tile breakdown based on tileClustersMax_
    // this will keep the current files open if the paths don't change
    bamClusterLoader_.open(bamFlowcellLayout_.getFlowcellId(), bamPath_);
    alignment::BclClusters::iterator clustersEnd = clusters_.cluster(0);
    clusters_.pf().clear();
    std::back_insert_iterator<std::vector<bool> > pfIt(clusters_.pf());
    const unsigned clustersLoaded = bamClusterLoader_.loadClusters(clustersToLoad,
        bamFlowcellLayout_.getReadNameLength(), bamFlowcellLayout_.getReadMetadataList(), clustersEnd, pfIt);
    ISAAC_THREAD_CERR<< "Loaded  " << clustersLoaded << " clusters of length " << clusterLength_ << std::endl;
    clusters_.reset(clusterLength_, clustersLoaded);
    return clustersLoaded;
}

void BamBaseCallsSource::discoverTiles(flowcell::TileMetadataList &ret)
{
    ret.clear();

    const unsigned clustersLoaded = loadNextTile();
    if (clustersLoaded)
    {
        const std::string &flowcellId = bamFlowcellLayout_.getFlowcellId();

        const flowcell::TileMetadata tileMetadata(
            flowcellId, bamFlowcellLayout_.getIndex(),
            currentTile_ + 1, 1,
            clustersLoaded,
            currentTile_);
        ++currentTile_;
        ret.push_back(tileMetadata);
        ISAAC_THREAD_CERR << "Generated bam tile: " <<
            oligo::bclToString(reinterpret_cast<const unsigned char *>(&*clusters_.cluster(0)), flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList())) << " " <<
            oligo::bclToString(reinterpret_cast<const unsigned char *>(&*clusters_.cluster(0)) +
                               flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList()) * (tileMetadata.getClusterCount() - 1),
                               flowcell::getTotalReadLength(bamFlowcellLayout_.getReadMetadataList())) << " " << tileMetadata << std::endl;
    }
}

void BamBaseCallsSource::loadClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "Loaded bam tile: " << tileMetadata << " with " << clusters_.getClusterCount() << " clusters" << std::endl;

    bclData.swap(clusters_);
}


void BackgroundBamBaseCallsSource::loadTilesThread()
{
    ISAAC_THREAD_CERR<< "BackgroundBamBaseCallsSource load thread started" << std::endl;
    std::unique_lock<std::mutex> lock(stateMutex_);
    while (!terminateRequested_)
    {
        unsigned loadedClustersCount = 0;
        ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::loadTilesThread 1" << std::endl;
        try
        {
            // make sure client thread can pick up previous tile while we're loading this one
            common::unlock_guard<std::unique_lock<std::mutex> > unlock(lock);
            loadedClustersCount = loadNextTile();
        }
        catch (...)
        {
            forceTermination_ = true;
            throw;
        }

        ISAAC_ASSERT_MSG(!loadedTile_, "Expected client to be ready for next tile");
        ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::loadTilesThread 2" << std::endl;

        if (loadedClustersCount)
        {
            loadedTile_ = ++loadingTile_;
            loadedClustersCount_ = loadedClustersCount;
            ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::loadTilesThread 4" << std::endl;
            stateChangeEvent_.notify_all();

            while (loadedTile_ && !terminateRequested_)
            {
                // last loaded tile has not been picked up
                stateChangeEvent_.wait(lock);
            }
        }
        else
        {
            loadedTile_ = 0;
            loadedClustersCount_ = 0;
            noMoreData_ = true;
            ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::loadTilesThread 5" << std::endl;
            stateChangeEvent_.notify_all();
            break;
        }
    }
    ISAAC_THREAD_CERR<< "BackgroundBamBaseCallsSource load thread terminated" << std::endl;
}


flowcell::TileMetadataList BackgroundBamBaseCallsSource::discoverTiles()
{
    ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::discoverTiles" << std::endl;
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
    if (loadedTile_)
    {
        ret.push_back(flowcell::TileMetadata(
            bamFlowcellLayout_.getFlowcellId(), bamFlowcellLayout_.getIndex(), loadedTile_,
            1, loadedClustersCount_, 0));
    }
    ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::discoverTiles done" << std::endl;
    return ret;
}

void BackgroundBamBaseCallsSource::loadClusters(
    const flowcell::TileMetadata &tileMetadata,
    alignment::BclClusters &bclData)
{
    ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::loadClusters" << std::endl;

    std::unique_lock<std::mutex> lock(stateMutex_);
    if (forceTermination_)
    {
        BOOST_THROW_EXCEPTION(common::ThreadingException("Terminating due to failures on other threads"));
    }

    ISAAC_ASSERT_MSG(tileMetadata.getFlowcellIndex() == bamFlowcellLayout_.getIndex(), "Unexpected tile requested " << tileMetadata);
    ISAAC_ASSERT_MSG(tileMetadata.getLane() == 1, "Unexpected tile lane requested: " << tileMetadata);
    ISAAC_ASSERT_MSG(tileMetadata.getTile() == loadedTile_, "Unexpected tile tile requested: " << tileMetadata);

    BamBaseCallsSource::loadClusters(tileMetadata, bclData);

    loadedTile_ = 0;
    stateChangeEvent_.notify_all();

    ISAAC_THREAD_CERR << "BackgroundBamBaseCallsSource::loadClusters done" << std::endl;
}


} // namespace alignWorkflow
} // namespace workflow
} // namespace isaac
