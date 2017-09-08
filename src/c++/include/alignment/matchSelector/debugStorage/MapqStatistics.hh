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
 ** \file DebugFragmentStorage.hh
 **
 ** \brief Compares alignment result with truth data in the read name and produces various statistics for accuracy debugging.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_MAPQ_STATISTICS_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_MAPQ_STATISTICS_HH

#include <atomic>

#include "alignment/matchSelector/FragmentStorage.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{
namespace debugStorage
{

class MapqStatistics
{
    const unsigned debugClassFilter_;
    const boost::filesystem::path outputMapqFilePath_;
    const boost::filesystem::path outputSmFilePath_;
    const boost::filesystem::path outputAsFilePath_;

    std::atomic<std::size_t> unalignedFragments_;
    std::atomic<std::size_t> downgradedMapq_;
    std::atomic<std::size_t> badDowngradedMapq_;
    std::atomic<std::size_t> downgradedSm_;
    std::atomic<std::size_t> badDowngradedSm_;
    std::atomic<std::size_t> downgradedAs_;
    std::atomic<std::size_t> badDowngradedAs_;
    std::array<std::atomic<std::size_t>, 2000> alignmentsBySm_;
    std::array<std::atomic<std::size_t>, 2000> badAlignmentsBySm_;
    std::array<std::atomic<std::size_t>, 2000> alignmentsByAs_;
    std::array<std::atomic<std::size_t>, 2000> badAlignmentsByAs_;
    std::array<std::atomic<std::size_t>, 256> alignmentsByMapq_;
    std::array<std::atomic<std::size_t>, 256> badAlignmentsByMapq_;
public:
    MapqStatistics(
        const unsigned debugClassFilter,
        const boost::filesystem::path &outputMapqFilePath,
        const boost::filesystem::path &outputSMFilePath,
        const boost::filesystem::path &outputAsFilePath);

    ~MapqStatistics();

    bool updateStat(
        const BamTemplate& bamTemplate,
        const std::size_t readIndex,
        const FragmentMetadata &fragment,
        const FragmentMetadata *mate);

private:
    bool updateMapqHistograms(
        const BamTemplate& bamTemplate,
        const std::size_t readIndex,
        const FragmentMetadata& fragment,
        const FragmentMetadata *mate);
};

} // namespace debugStorage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_STORAGE_MAPQ_STATISTICS_HH
