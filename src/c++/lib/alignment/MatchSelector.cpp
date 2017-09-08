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
 ** \file MatchSelector.cpp
 **
 ** Component to select the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#include <numeric>
#include <fstream>
#include <cerrno>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/thread.hpp>

#include "alignment/Mismatch.hh"
#include "alignment/MatchSelector.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "common/FastIo.hh"
#include "reference/Contig.hh"
#include "reference/ContigLoader.hh"

#include "alignment/matchSelector/MatchSelectorStatsXml.hh"

namespace isaac
{
namespace alignment
{


std::vector<SequencingAdapterList> generateSequencingAdapters(const flowcell::BarcodeMetadataList &barcodeMetadataList)
{
    std::vector<SequencingAdapterList> ret(barcodeMetadataList.size());

    std::vector<SequencingAdapterList>::iterator barcodeAdaptersIterator = ret.begin();
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcode, barcodeMetadataList)
    {
        BOOST_FOREACH(const flowcell::SequencingAdapterMetadata &adapter, barcode.getAdapters())
        {
            barcodeAdaptersIterator->push_back(SequencingAdapter(adapter));
        }
        ++barcodeAdaptersIterator;
    }

    return ret;
}

MatchSelector::MatchSelector(
        const unsigned int maxThreadCount,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const reference::NumaContigLists &contigLists,
        const std::size_t candidateMatchesMax,
        const unsigned repeatThreshold,
        const std::vector<std::size_t> &clusterIdList,
        const int mateDriftRange,
        const TemplateLengthStatistics &userTemplateLengthStatistics,
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
        const unsigned detectTemplateBlockSize
    )
    : computeThreads_(maxThreadCount),
      tileMetadataList_(),//(tileMetadataList),
      barcodeMetadataList_(barcodeMetadataList),
      flowcellLayoutList_(flowcellLayoutList),
      contigLists_(contigLists),
      repeatThreshold_(repeatThreshold),
      clusterIdList_(clusterIdList),
      mapqThreshold_(mapqThreshold),
      pfOnly_(pfOnly),
      collectCycleStats_(collectCycleStats),
      baseQualityCutoff_(baseQualityCutoff),
      keepUnaligned_(keepUnaligned),
      clipSemialigned_(clipSemialigned),
      clipOverlapping_(clipOverlapping),
      barcodeSequencingAdapters_(generateSequencingAdapters(barcodeMetadataList_)),
      allStats_(),//(tileMetadataList_.size(), matchSelector::MatchSelectorStats(barcodeMetadataList_)),
      threadStats_(computeThreads_.size(), matchSelector::MatchSelectorStats(collectCycleStats_, barcodeMetadataList_)),
      threadCluster_(computeThreads_.size(),
                     Cluster(flowcell::getMaxReadLength(flowcellLayoutList_) +
                             flowcell::getMaxBarcodeLength(flowcellLayoutList_))),
      threadTemplateBuilders_(computeThreads_.size()),
      threadSemialignedEndsClippers_(clipSemialigned_ ? computeThreads_.size() : 0),
      threadOverlappingEndsClippers_(computeThreads_.size()),
      restOfGenomeCorrections_(barcodeMetadataList_.size()),
      templateDetector_(
          computeThreads_,
          barcodeMetadataList_,
          flowcellLayoutList_,
          contigLists_,
          threadTemplateBuilders_,
          mateDriftRange,
          userTemplateLengthStatistics,
          perTileTls,
          detectTemplateBlockSize)
{
    ISAAC_TRACE_STAT("Constructing match selector");
    while(threadTemplateBuilders_.size() < computeThreads_.size())
    {
        threadTemplateBuilders_.push_back(new TemplateBuilder(collectCycleStats_,
                                                              flowcellLayoutList_,
                                                              candidateMatchesMax,
                                                              repeatThreshold_,
                                                              seedLength,
                                                              maxSeedsPerMatch,
                                                              matchFinderTooManyRepeats,
                                                              matchFinderWayTooManyRepeats,
                                                              matchFinderShadowSplitRepeats,
                                                              scatterRepeats,
                                                              rescueShadows,
                                                              trimPEAdapters,
                                                              anchorMate,
                                                              gappedMismatchesMax,
                                                              smitWatermanGapsMax,
                                                              smartSmithWaterman,
                                                              smithWatermanGapSizeMax,
                                                              splitAlignments,
                                                              alignmentCfg,
                                                              dodgyAlignmentScore, anomalousPairHandicap, reserveBuffers));
    }
    ISAAC_TRACE_STAT("Constructed match selector");
}

void MatchSelector::dumpStats(const boost::filesystem::path &statsXmlPath)
{
    std::for_each(allStats_.begin(), allStats_.end(), boost::bind(&matchSelector::MatchSelectorStats::finalize, _1));

    std::ofstream os(statsXmlPath.string().c_str());
    if (!os) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "ERROR: Unable to open file for writing: " + statsXmlPath.string()));
    }

    matchSelector::MatchSelectorStatsXml statsXml(
        collectCycleStats_, flowcellLayoutList_, barcodeMetadataList_, tileMetadataList_, allStats_);
    statsXml.serialize(os);
}

/**
 * \return false if cluster is unaligned
 */
template <typename MatchFinderT>
matchSelector::TemplateAlignmentType MatchSelector::alignCluster(
    const reference::ContigList& barcodeContigList,
    const flowcell::ReadMetadataList& tileReads,
    const SequencingAdapterList& sequencingAdapters,
    const TemplateLengthStatistics& templateLengthStatistics,
    const uint64_t barcodeIndex,
    const MatchFinderT &matchFinder,
    const RestOfGenomeCorrection& restOfGenomeCorrection,
    const unsigned threadNumber, TemplateBuilder& ourThreadTemplateBuilder,
    const Cluster& cluster, BamTemplate& bamTemplate,
    matchSelector::MatchSelectorStats& stats,
    matchSelector::FragmentStorage &fragmentStorage)
{
    const templateBuilder::AlignmentType res = ourThreadTemplateBuilder.buildTemplate(
        barcodeContigList, restOfGenomeCorrection, tileReads,
        sequencingAdapters, cluster, templateLengthStatistics, true, matchFinder, bamTemplate);
    // build the fragments for that cluster
    if (templateBuilder::Normal == res)
    {
        if (clipSemialigned_)
        {
            threadSemialignedEndsClippers_[threadNumber].reset();
            threadSemialignedEndsClippers_[threadNumber].clip(barcodeContigList, bamTemplate);
        }
        if (clipOverlapping_)
        {
            threadOverlappingEndsClippers_[threadNumber].reset();
            threadOverlappingEndsClippers_[threadNumber].clip(barcodeContigList, bamTemplate);
        }
    }

    if (templateBuilder::Normal == res || keepUnaligned_)
    {
        if (bamTemplate.filterLowQualityFragments(mapqThreshold_))
        {
            return matchSelector::Qc;
        }

        fragmentStorage.store(bamTemplate, barcodeIndex, threadNumber);
    }
    return templateBuilder::Nm == res ? matchSelector::NmNm : templateBuilder::Rm == res ? matchSelector::Rm : matchSelector::Qc;
}

template <typename MatchFinderT>
void MatchSelector::alignThread(
    const unsigned threadNumber,
    const flowcell::TileMetadata & tileMetadata,
    const matchFinder::ClusterInfos &clusterInfos,
    unsigned &threadClusterId,
    const MatchFinderT &matchFinder,
    const BclClusters &bclData,
    const std::vector<TemplateLengthStatistics> & templateLengthStatistics,
    matchSelector::FragmentStorage &fragmentStorage)
{
    Cluster &ourThreadCluster = threadCluster_[threadNumber];
    TemplateBuilder &ourThreadTemplateBuilder = threadTemplateBuilders_.at(threadNumber);
    matchSelector::MatchSelectorStats &ourThreadStats = threadStats_.at(threadNumber);

    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const flowcell::ReadMetadataList &tileReads = flowcell.getReadMetadataList();
    const std::size_t barcodeLength = flowcell.getBarcodeLength();
    const unsigned readNameLength = flowcell.getReadNameLength();

    const reference::ContigLists &threadContigLists = contigLists_.threadNodeContainer();

    boost::unique_lock<boost::mutex> lock(mutex_);

    while (tileMetadata.getClusterCount() != threadClusterId)
    {
        const unsigned clustersBegin = threadClusterId;
        static const unsigned clustersAtATime = CLUSTERS_AT_A_TIME;
        threadClusterId += std::min(clustersAtATime, tileMetadata.getClusterCount() - threadClusterId);
        const unsigned clustersEnd = threadClusterId;
        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            for (unsigned clusterId = clustersBegin; clustersEnd != clusterId; ++clusterId)
            {
                if (!clusterIdList_.empty() && clusterIdList_.end() == std::find(clusterIdList_.begin(), clusterIdList_.end(), clusterId))
                {
                    continue;
                }
                const flowcell::BarcodeMetadata &barcodeMetadata = barcodeMetadataList_[clusterInfos[clusterId].getBarcodeIndex()];

                // uninitialize cluster in case it does not get stored in as storage that buffers data
                // not relevant anymore as BufferingFragmentStorage is gone
                fragmentStorage.reset(clusterId, 2 == tileReads.size());

                // initialize the cluster with the bcl data
                ourThreadCluster.init(tileReads, bclData.cluster(clusterId), tileMetadata.getIndex(), clusterId,
                                      bclData.xy(clusterId), bclData.pf(clusterId), barcodeLength, readNameLength);
                BamTemplate bamTemplate(tileReads, ourThreadCluster);

                matchSelector::TemplateAlignmentType result = matchSelector::Filtered;
                if (!barcodeMetadata.isUnmappedReference())
                {
                    const reference::ContigList &barcodeContigList = threadContigLists.at(barcodeMetadata.getReferenceIndex());
                    const SequencingAdapterList &sequencingAdapters = barcodeSequencingAdapters_.at(barcodeMetadata.getIndex());

                    ISAAC_ASSERT_MSG(clusterId < tileMetadata.getClusterCount(), "Cluster ids are expected to be 0-based within the tile.");

                    trimLowQualityEnds(ourThreadCluster, baseQualityCutoff_);

                    // if pfOnly_ is set, this non-pf cluster will not be reported as a regularly-processed one.
                    // if match list begins with noMatchReferencePosition, then this cluster does not have any matches at all. This is
                    // because noMatchReferencePosition has the highest possible contig number and sort will put it to the end of match list
                    // In either case report it as skipped to ensure statistics consistency
                    if (!pfOnly_ || bclData.pf(clusterId))
                    {
                        result = alignCluster(
                            barcodeContigList, tileReads, sequencingAdapters,
                            templateLengthStatistics[barcodeMetadata.getIndex()], barcodeMetadata.getIndex(), matchFinder,
                            restOfGenomeCorrections_[barcodeMetadata.getIndex()],
                            threadNumber, ourThreadTemplateBuilder, ourThreadCluster, bamTemplate, ourThreadStats,
                            fragmentStorage);
                    }
                }
                ourThreadStats.recordTemplate(
                    tileReads, templateLengthStatistics[barcodeMetadata.getIndex()],
                    bamTemplate, barcodeMetadata.getIndex(), result);
            }
        }
    }
}

template <typename MatchFinderT>
void MatchSelector::parallelSelect(
    alignment::matchFinder::TileClusterInfo &tileClusterInfo,
    std::vector<TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
    const flowcell::TileMetadata &tileMetadata,
    const MatchFinderT &matchFinder,
    const BclClusters &bclData,
    matchSelector::FragmentStorage &fragmentStorage)
{
    const reference::ContigLists &threadContigLists = contigLists_.threadNodeContainer();
    std::for_each(threadStats_.begin(), threadStats_.end(), boost::bind(&matchSelector::MatchSelectorStats::reset, _1));

    ISAAC_THREAD_CERR << "Resizing fragment storage for " <<  tileMetadata.getClusterCount() << " clusters " << std::endl;
    fragmentStorage.resize(tileMetadata.getClusterCount());
    ISAAC_THREAD_CERR << "Resizing fragment storage done for " <<  tileMetadata.getClusterCount() << " clusters " << std::endl;

    // Recompute relevant genome corrections
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const flowcell::ReadMetadataList &tileReads = flowcell.getReadMetadataList();
    BOOST_FOREACH(const flowcell::BarcodeMetadata &barcodeMetadata, barcodeMetadataList_)
    {
        if (tileMetadata.getLane() == barcodeMetadata.getLane() && !barcodeMetadata.isUnmappedReference())
        {
            const reference::ContigList &barcodeContigList = threadContigLists.at(barcodeMetadata.getReferenceIndex());
            restOfGenomeCorrections_[barcodeMetadata.getIndex()] = RestOfGenomeCorrection(barcodeContigList, tileReads);
        }
    }

    templateDetector_.determineTemplateLengths(
        tileMetadata, tileClusterInfo.at(tileMetadata.getIndex()),
        restOfGenomeCorrections_,
        matchFinder, bclData, barcodeTemplateLengthStatistics, threadStats_[0]);

    ISAAC_THREAD_CERR << "Selecting matches on " <<  computeThreads_.size() << " threads for " << tileMetadata << "\n" << std::endl;
    unsigned clusterId = 0;
    computeThreads_.execute(boost::bind(&MatchSelector::alignThread<MatchFinderT>, this, _1,
                                        boost::ref(tileMetadata),
                                        boost::ref(tileClusterInfo.at(tileMetadata.getIndex())),
                                        boost::ref(clusterId),
                                        boost::ref(matchFinder),
                                        boost::ref(bclData),
                                        boost::cref(barcodeTemplateLengthStatistics),
                                        boost::ref(fragmentStorage)));

    ISAAC_THREAD_CERR << "Selecting matches done on " <<  computeThreads_.size() << " threads for " << clusterId << " clusters of " << tileMetadata  << std::endl;

    BOOST_FOREACH(const matchSelector::MatchSelectorStats &threadStats, threadStats_)
    {
        allStats_.at(tileMetadata.getIndex()) += threadStats;
    }
}

void MatchSelector::reserveMemory(
    const flowcell::TileMetadataList &tileMetadataList)
{
    for (const flowcell::TileMetadata &tileMetadata : tileMetadataList)
    {
        tileMetadataList_.resize(std::max<std::size_t>(tileMetadata.getIndex() + 1, tileMetadataList_.size()));
        tileMetadataList_.at(tileMetadata.getIndex()) = tileMetadata;
        allStats_.resize(std::max<std::size_t>(tileMetadata.getIndex() + 1, allStats_.size()), matchSelector::MatchSelectorStats(collectCycleStats_, barcodeMetadataList_));
    }
}

template <typename KmerT> struct InstantiateTemplates : MatchSelector
{
    typedef ClusterHashMatchFinder<reference::ReferenceHash<KmerT, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > > MatchFinderT;
    void parallelSelectInstance(alignment::matchFinder::TileClusterInfo &tileClusterInfo,
                                std::vector<TemplateLengthStatistics> &barcodeTemplateLengthStatistics,
                                const flowcell::TileMetadata &tileMetadata,
                                const MatchFinderT &matchFinder,
                                const BclClusters &bclData,
                                matchSelector::FragmentStorage &fragmentStorage)
    {
        MatchSelector::parallelSelect(tileClusterInfo, barcodeTemplateLengthStatistics, tileMetadata, matchFinder, bclData, fragmentStorage);
    }
};

template struct InstantiateTemplates<oligo::BasicKmerType<10> >;
template struct InstantiateTemplates<oligo::BasicKmerType<11> >;
template struct InstantiateTemplates<oligo::BasicKmerType<12> >;
template struct InstantiateTemplates<oligo::BasicKmerType<13> >;
template struct InstantiateTemplates<oligo::BasicKmerType<14> >;
template struct InstantiateTemplates<oligo::BasicKmerType<15> >;
template struct InstantiateTemplates<oligo::BasicKmerType<16> >;
template struct InstantiateTemplates<oligo::BasicKmerType<17> >;
template struct InstantiateTemplates<oligo::BasicKmerType<18> >;
template struct InstantiateTemplates<oligo::BasicKmerType<19> >;
template struct InstantiateTemplates<oligo::BasicKmerType<20> >;

} // namespace alignemnt
} // namespace isaac
