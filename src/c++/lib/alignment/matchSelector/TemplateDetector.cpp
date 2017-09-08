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
 ** \file TemplateDetector.cpp
 **
 ** Component to select the best matches among all possible candidates.
 **
 ** \author Come Raczy
 **/

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "alignment/templateBuilder/FragmentSequencingAdapterClipper.hh"
#include "alignment/matchSelector/TemplateDetector.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{


TemplateDetector::TemplateDetector(
    common::UnsafeThreadVector &computeThreads,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const reference::NumaContigLists &contigLists,
    boost::ptr_vector<TemplateBuilder> &threadTemplateBuilders,
    const int mateDriftRange,
    const TemplateLengthStatistics &userTemplateLengthStatistics,
    const bool perTileTls,
    const unsigned detectTemplateBlockSize)
: computeThreads_(computeThreads),
  barcodeMetadataList_(barcodeMetadataList),
  flowcellLayoutList_(flowcellLayoutList),
  contigLists_(contigLists),
  userTemplateLengthStatistics_(userTemplateLengthStatistics),
  perTileTls_(perTileTls),
  threadCluster_(computeThreads_.size(),
                 Cluster(flowcell::getMaxReadLength(flowcellLayoutList_) +
                         flowcell::getMaxBarcodeLength(flowcellLayoutList_))),
  threadTemplateBuilders_(threadTemplateBuilders),
  templateLengthDistributions_(barcodeMetadataList_.size(), TemplateLengthDistribution(detectTemplateBlockSize, mateDriftRange)),
  unprocessedClusterId_(0),
  pendingClusterId_(0),
  detectTemplateBlockSize_(detectTemplateBlockSize)
{
}

template <typename MatchFinderT>
void TemplateDetector::collectModels(
    const unsigned clusterRangeBegin,
    const unsigned clusterRangeEnd,
    const BclClusters& bclData,
    const matchFinder::ClusterInfos& clusterInfos,
    const flowcell::ReadMetadataList& tileReads,
    const flowcell::TileMetadata& tileMetadata,
    const unsigned barcodeLength,
    const unsigned readNameLength,
    std::size_t& statsToBuild,
    std::vector<alignment::TemplateLengthStatistics>& templateLengthStatistics,
    Cluster& ourThreadCluster,
    TemplateBuilder& templateBuilder,
    const MatchFinderT& matchFinder,
    common::StaticVector<BarcodeAlignmentModel, CLUSTERS_AT_A_TIME>& threadBarcodeModels)
{
    const reference::ContigLists &threadContigLists = contigLists_.threadNodeContainer();
    for (unsigned clusterId = clusterRangeBegin; statsToBuild && clusterId != clusterRangeEnd; ++clusterId)
    {
        if (!bclData.pf(clusterId))
        {
            continue;
        }
        const unsigned barcodeIndex = clusterInfos[clusterId].getBarcodeIndex();
        if (barcodeMetadataList_[barcodeIndex].isUnmappedReference() || templateLengthStatistics[barcodeIndex].isStable())
        {
            continue;
        }
        const unsigned barcodeReference = barcodeMetadataList_[barcodeIndex].getReferenceIndex();

        ourThreadCluster.init(
            tileReads, bclData.cluster(clusterId), tileMetadata.getIndex(), clusterId,
            bclData.xy(clusterId), true, barcodeLength, readNameLength);

        SequencingAdapterList emptyList;
        templateBuilder::FragmentSequencingAdapterClipper adapterClipper(emptyList);
        const templateBuilder::AlignmentType alignmentType = templateBuilder.buildFragments(
                threadContigLists.at(barcodeReference),
                tileReads, adapterClipper,
                matchFinder, MATCH_FINDER_TOO_MANY_REPEATS, ourThreadCluster, false);

        if (templateBuilder::Normal == alignmentType)
        {
            const TemplateLengthDistribution::AlignmentModel am = TemplateLengthDistribution::getAlignmentModel(
                templateBuilder.getFragments()[0], templateBuilder.getFragments()[1]);
            if (TemplateLengthStatistics::InvalidAlignmentModel != am.first)
            {
                threadBarcodeModels.push_back(std::make_pair(barcodeIndex, am));
            }
        }
    }
}

template <typename MatchFinderT>
void TemplateDetector::templateLengthThread(
    const unsigned threadNumber,
    const flowcell::TileMetadata& tileMetadata, const BclClusters& bclData,
    const matchFinder::ClusterInfos& clusterInfos,
    const MatchFinderT& matchFinder,
    std::size_t &statsToBuild,
    std::vector<alignment::TemplateLengthStatistics>& templateLengthStatistics)
{
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const flowcell::ReadMetadataList &tileReads = flowcell.getReadMetadataList();

    const unsigned barcodeLength = flowcell.getBarcodeLength();
    const unsigned readNameLength = flowcell.getReadNameLength();
    TemplateBuilder& ourThreadTemplateBuilder = threadTemplateBuilders_.at(threadNumber);
    Cluster& ourThreadCluster = threadCluster_.at(threadNumber);
    const alignment::BclClusterFields<> fieldsParser(tileReads, flowcell.getBarcodeLength());


    boost::unique_lock<boost::mutex> lock(mutex_);
    while (statsToBuild && tileMetadata.getClusterCount() != unprocessedClusterId_)
    {
        common::StaticVector<BarcodeAlignmentModel, CLUSTERS_AT_A_TIME> threadBarcodeModels;
        const unsigned clusterRangeBegin = unprocessedClusterId_;
        static const unsigned clustersAtATime = CLUSTERS_AT_A_TIME;
        const unsigned clusterRangeEnd = clusterRangeBegin + std::min(clustersAtATime, tileMetadata.getClusterCount() - clusterRangeBegin);
        unprocessedClusterId_ = clusterRangeEnd;

        {
            common::unlock_guard<boost::unique_lock<boost::mutex> > unlock(lock);
            collectModels(clusterRangeBegin, clusterRangeEnd, bclData,
                          clusterInfos, tileReads, tileMetadata, barcodeLength,
                          readNameLength, statsToBuild,
                          templateLengthStatistics, ourThreadCluster,
                          ourThreadTemplateBuilder, matchFinder,
                          threadBarcodeModels);
        }

        // now that the lock has been reacquired, we need to wait for our turn to update the distributions
        while (statsToBuild && pendingClusterId_ != clusterRangeBegin)
        {
            stateChangedCondition_.wait(lock);
        }

        for(const BarcodeAlignmentModel &bm : threadBarcodeModels)
        {
            const unsigned barcodeIndex = bm.first;
            if (!templateLengthDistributions_[barcodeIndex].isStable())
            {
                templateLengthDistributions_[barcodeIndex].appendModel(bm.second.first, bm.second.second);
                if (templateLengthDistributions_[barcodeIndex].isStable())
                {
                    templateLengthStatistics[barcodeIndex] = templateLengthDistributions_[barcodeIndex].getStatistics();

                    // just stabilized another one which was not stable. Terminate early when there is nothing left to build
                    --statsToBuild;
                    ISAAC_THREAD_CERR << "Determining template length done for " << tileMetadata << ", " << barcodeMetadataList_[barcodeIndex] << ":" << templateLengthStatistics[barcodeIndex] << " " << statsToBuild << " to go" << std::endl;
                }
            }
        }

        pendingClusterId_ = clusterRangeEnd;
        stateChangedCondition_.notify_all();
    }
}

template <typename MatchFinderT>
void TemplateDetector::determineTemplateLengths(
    const flowcell::TileMetadata &tileMetadata,
    const matchFinder::ClusterInfos &clusterInfos,
    const std::vector<RestOfGenomeCorrection> &restOfGenomeCorrections,
    const MatchFinderT &matchFinder,
    const BclClusters &bclData,
    std::vector<alignment::TemplateLengthStatistics> &templateLengthStatistics,
    matchSelector::MatchSelectorStats &stats)
{
    const flowcell::Layout &flowcell = flowcellLayoutList_.at(tileMetadata.getFlowcellIndex());
    const flowcell::ReadMetadataList &tileReads = flowcell.getReadMetadataList();

    ISAAC_ASSERT_MSG(2 >= tileReads.size(), "only single-ended and paired reads are supported");

    if (2 != tileReads.size())
    {
        ISAAC_THREAD_CERR << "Using unstable template-length statistics for single-ended data" << std::endl;
        return;
    }

    if (userTemplateLengthStatistics_.isStable())
    {
        ISAAC_THREAD_CERR << "Using user-defined template-length statistics: " << userTemplateLengthStatistics_ << std::endl;
        std::fill(templateLengthStatistics.begin(), templateLengthStatistics.end(), userTemplateLengthStatistics_);
        return;
    }

    if (perTileTls_)
    {
        ISAAC_THREAD_CERR << "Resetting template length for " << tileMetadata << std::endl;
        std::for_each(
            templateLengthStatistics.begin(), templateLengthStatistics.end(),
            [&](TemplateLengthStatistics &btls)
            {
                const flowcell::BarcodeMetadata &barcodeMetadata = barcodeMetadataList_[std::distance(&templateLengthStatistics.front(), &btls)];
                if (barcodeMetadata.getLane() == tileMetadata.getLane() &&
                    barcodeMetadata.getFlowcellIndex() == flowcell.getIndex())
                {
                    btls.clear();
                }
            });
    }

    std::for_each(templateLengthDistributions_.begin(), templateLengthDistributions_.end(), boost::bind(&TemplateLengthDistribution::clear, _1));

    std::size_t statsToBuild = 0;
    std::size_t barcode = 0;
    BOOST_FOREACH(TemplateLengthStatistics &tls, templateLengthStatistics)
    {
        statsToBuild +=
            !tls.isStable() &&
            barcodeMetadataList_[barcode].getLane() == tileMetadata.getLane() &&
            barcodeMetadataList_[barcode].getFlowcellIndex() == flowcell.getIndex();
        ++barcode;
    }

    ISAAC_THREAD_CERR << "Determining template length statistics for " << statsToBuild << " barcodes on " << tileMetadata << std::endl;

    if (statsToBuild)
    {
        unprocessedClusterId_ = 0;
        pendingClusterId_ = 0;

        computeThreads_.execute(boost::bind(&TemplateDetector::templateLengthThread<MatchFinderT>, this, _1,
                                            boost::ref(tileMetadata),
                                            boost::ref(bclData),
                                            boost::ref(clusterInfos),
                                            boost::ref(matchFinder),
                                            boost::ref(statsToBuild),
                                            boost::ref(templateLengthStatistics)));
    }

    barcode = 0;
    BOOST_FOREACH(TemplateLengthDistribution &templateLengthDistribution, templateLengthDistributions_)
    {
        const flowcell::BarcodeMetadata &barcodeMetadata = barcodeMetadataList_[barcode];

        if (barcodeMetadata.getLane() == tileMetadata.getLane() &&
            barcodeMetadata.getFlowcellIndex() == flowcell.getIndex())
        {
            if (!templateLengthStatistics[barcode].isStable() && !templateLengthDistribution.isStable())
            {
                templateLengthDistribution.finalize();
                templateLengthStatistics[barcode] = templateLengthDistribution.getStatistics();
                ISAAC_THREAD_CERR << "Unstable template length finalized for " << tileMetadata << ", " << barcodeMetadata << ":" << templateLengthStatistics[barcode] << std::endl;
            }
            stats.recordTemplateLengthStatistics(barcodeMetadata, templateLengthStatistics[barcode]);
        }
        ++barcode;
    }
}

template <typename KmerT> struct InstantiateTemplates : TemplateDetector
{
    typedef ClusterHashMatchFinder<reference::ReferenceHash<KmerT, common::NumaAllocator<void, common::numa::defaultNodeInterleave> >> MatchFinderT;

    void determineTemplateLengths(
        const flowcell::TileMetadata &tileMetadata,
        const matchFinder::ClusterInfos &clusterInfos,
        const std::vector<RestOfGenomeCorrection> &restOfGenomeCorrection,
        const MatchFinderT &matchFinder,
        const BclClusters &bclData,
        std::vector<alignment::TemplateLengthStatistics> &templateLengthStatistics,
        matchSelector::MatchSelectorStats &stats)
    {
        TemplateDetector::determineTemplateLengths(tileMetadata, clusterInfos, restOfGenomeCorrection, matchFinder, bclData,
                templateLengthStatistics, stats);
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

} // namespace templateDetector
} // namespace alignemnt
} // namespace isaac
