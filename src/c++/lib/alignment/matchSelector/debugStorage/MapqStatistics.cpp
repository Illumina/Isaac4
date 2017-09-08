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
 ** \file DebugFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>

#include "alignment/matchSelector/debugStorage/MapqStatistics.hh"
#include "alignment/matchSelector/debugStorage/Utils.hh"
#include "common/Debug.hh"
#include "common/FastIo.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{
namespace debugStorage
{

//std::array<std::array<std::atomic<std::size_t>, 6 >, 6> badAlignmentsByGoodAndBadSeedCount;
//std::array<std::array<std::atomic<std::size_t>, 6 >, 6> goodAlignmentsByGoodAndBadSeedCount;

//void traceMatrix(
//    const std::array<std::array<std::atomic<std::size_t>, 6>, 6>& alignmentsByGoodAndBadSeedCount,
//    std::string matrix)
//{
//    matrix += "\n";
//    for (const auto& line : alignmentsByGoodAndBadSeedCount)
//    {
//        for (const auto& cell : line)
//        {
//            matrix += "\t";
//            matrix += std::to_string(cell);
//        }
//        matrix += "\n";
//    }
//    std::cout << matrix << std::endl;
//}

template <typename AlignmentsByScore>
void traceQQ(
    std::ostream &os,
    const std::atomic<std::size_t> &unaligned,
    const std::atomic<std::size_t> &downgraded,
    const std::atomic<std::size_t> &downgradedBad,
    const AlignmentsByScore &alignmentsByScore,
    const AlignmentsByScore &badAlignmentsByScore)
{
    os << "-2\t0\t0\t" << unaligned << std::endl;
    if (downgraded)
    {
        os << "-1\t" << -10 * log10(double(downgradedBad) / double(downgraded)) << "\t" << downgradedBad << "\t" << downgraded << std::endl;
    }
    const std::size_t toTrace =
        std::distance(
            alignmentsByScore.begin(),
            std::find_if(alignmentsByScore.rbegin(), alignmentsByScore.rend(),
                [](const std::size_t c){return c;}).base());
    for (std::size_t i= 0; i <= toTrace; ++i)
    {
        if (alignmentsByScore[i])
        {
            os << i << "\t" << -10 * log10(double(badAlignmentsByScore[i]) / double(alignmentsByScore[i])) << "\t" << badAlignmentsByScore[i] << "\t" << alignmentsByScore[i] << std::endl;
        }
    }
}

MapqStatistics::~MapqStatistics()
{
//    traceMatrix(badAlignmentsByGoodAndBadSeedCount, "badAlignmentsByGoodAndBadSeedCount");
//    traceMatrix(goodAlignmentsByGoodAndBadSeedCount, "goodAlignmentsByGoodAndBadSeedCount");
//
    if (alignmentsBySm_.end() != std::find_if(alignmentsBySm_.begin(), alignmentsBySm_.end(), [](const std::atomic<std::size_t> &c){return 0 != c;}) ||
        alignmentsByAs_.end() != std::find_if(alignmentsByAs_.begin(), alignmentsByAs_.end(), [](const std::atomic<std::size_t> &c){return 0 != c;}) ||
        alignmentsByMapq_.end() != std::find_if(alignmentsByMapq_.begin(), alignmentsByMapq_.end(), [](const std::atomic<std::size_t> &c){return 0 != c;}))
    {
        {
            ISAAC_THREAD_CERR << "Storing SM QQ stats in " << outputSmFilePath_.c_str() << std::endl;
            std::ofstream smTsv(outputSmFilePath_.c_str());
            traceQQ(smTsv, unalignedFragments_, downgradedSm_, badDowngradedSm_, alignmentsBySm_, badAlignmentsBySm_);
        }

        {
            ISAAC_THREAD_CERR << "Storing AS QQ stats in " << outputAsFilePath_.c_str() << std::endl;
            std::ofstream asTsv(outputAsFilePath_.c_str());
            traceQQ(asTsv, unalignedFragments_, downgradedAs_, badDowngradedAs_, alignmentsByAs_, badAlignmentsByAs_);
        }

        {
            ISAAC_THREAD_CERR << "Storing MAPQ QQ stats in " << outputMapqFilePath_.c_str() << std::endl;
            std::ofstream mapqTsv(outputMapqFilePath_.c_str());
            traceQQ(mapqTsv, unalignedFragments_, downgradedMapq_, badDowngradedMapq_, alignmentsByMapq_, badAlignmentsByMapq_);
        }
    }
}
//
//std::size_t getGoodSeedCount(const FragmentMetadata &fragment)
//{
//    return fragment.isAligned();
//}
//
//std::size_t getBadSeedCount(const FragmentMetadata &fragment)
//{
//    return !fragment.isAligned();
//}

bool MapqStatistics::updateMapqHistograms(
    const BamTemplate& bamTemplate,
    const std::size_t readIndex,
    const FragmentMetadata& fragment,
    const FragmentMetadata *mate)
{
    bool ret = true;

    const bool correct = alignsCorrectly(readIndex, fragment);

    if (fragment.hasAlignmentScore())
    {
        ++alignmentsBySm_[std::min<unsigned>(fragment.getAlignmentScore(), alignmentsBySm_.size() - 1)];
//        ret = true;
    }
    else
    {
        ++downgradedSm_;
        badDowngradedSm_ += correct;
    }

    if (bamTemplate.hasAlignmentScore() && fragment.hasAlignmentScore() && mate->hasAlignmentScore() && bamTemplate.isProperPair())
    {
        ++alignmentsByAs_[std::min<unsigned>(bamTemplate.getAlignmentScore(), alignmentsByAs_.size() - 1)];
    }
    else
    {
        ++downgradedAs_;
        badDowngradedAs_ += correct;
    }

    if (fragment.hasMapQ())
    {
        ++alignmentsByMapq_[std::min<unsigned>(fragment.mapQ, alignmentsByMapq_.size() - 1)];
    }
    else
    {
        ++downgradedMapq_;
        badDowngradedMapq_ += correct;
    }
//    ++alignmentsByMapq_[std::min<unsigned>(mapQ, alignmentsByMapq_.size() - 1)];
    if (!correct)
    {
        if (fragment.hasAlignmentScore())
        {
            ++badAlignmentsBySm_[std::min<unsigned >(fragment.getAlignmentScore(), badAlignmentsBySm_.size() - 1)];
        }
//        if (bamTemplate.hasAlignmentScore())
        if (bamTemplate.hasAlignmentScore() && fragment.hasAlignmentScore() && mate->hasAlignmentScore() && bamTemplate.isProperPair())
        {
            ++badAlignmentsByAs_[std::min<unsigned >(bamTemplate.getAlignmentScore(), badAlignmentsByAs_.size() - 1)];
        }
        if (fragment.hasMapQ())
        {
            ++badAlignmentsByMapq_[std::min<unsigned >(fragment.mapQ, badAlignmentsByMapq_.size() - 1)];
        }
    }
    return ret;
}

bool MapqStatistics::updateStat(
    const BamTemplate& bamTemplate,
    const std::size_t readNumber,
    const FragmentMetadata &fragment,
    const FragmentMetadata *mate)
{
    if (/*debugClassFilter_ != bamTemplate.debugClass_||*/ !fragment.isAligned())
    {
        ++unalignedFragments_;
    }
    else
    {
//        if (fragment.mismatchCount || 2 != fragment.secondBestMismatchDelta)//fragment.dodgy)
        return updateMapqHistograms(bamTemplate, readNumber, fragment, mate);
    }

    return false;
}


MapqStatistics::MapqStatistics(
    const unsigned debugClassFilter,
    const boost::filesystem::path &outputMapqFilePath,
    const boost::filesystem::path &outputSmFilePath,
    const boost::filesystem::path &outputAsFilePath):
    debugClassFilter_(debugClassFilter),
    outputMapqFilePath_(outputMapqFilePath),
    outputSmFilePath_(outputSmFilePath),
    outputAsFilePath_(outputAsFilePath),
    unalignedFragments_(0),
    downgradedMapq_(0),
    badDowngradedMapq_(0),
    downgradedSm_(0),
    badDowngradedSm_(0),
    downgradedAs_(0),
    badDowngradedAs_(0)
{
    for (auto &a : alignmentsBySm_) {a = 0;}
    for (auto &a : badAlignmentsBySm_) {a = 0;}
    for (auto &a : alignmentsByAs_) {a = 0;}
    for (auto &a : badAlignmentsByAs_) {a = 0;}
    for (auto &a : alignmentsByMapq_) {a = 0;}
    for (auto &a : badAlignmentsByMapq_) {a = 0;}
}

} // namespace debugStroage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac
