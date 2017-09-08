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
 ** \file FragmentMetadataTileStatsAdapter.hh
 **
 ** \brief Conversion from FragmentMetadata to the interface suitable for MatchSelectorTileStats::recordFragment
 ** generation.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_METADATA_TILE_STATS_ADAPTER_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_METADATA_TILE_STATS_ADAPTER_HH

#include <numeric>

#include <boost/noncopyable.hpp>

#include "alignment/FragmentMetadata.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class FragmentMetadataTileStatsAdapter: boost::noncopyable
{
    const FragmentMetadata &fragment_;
public:
    FragmentMetadataTileStatsAdapter(const FragmentMetadata &fragment):
        fragment_(fragment){}

    uint64_t getYield() const {return fragment_.getReadLength();}

    uint64_t getYieldQ30() const
    {
        const Read &read = fragment_.getRead();
        return std::count_if( read.getForwardQuality().begin(), read.getForwardQuality().end(),
                              [](unsigned char uc){return uc >= 30u;});
//                              boost::bind(&boost::cref<unsigned char>, _1) >= 30u);
    }

    uint64_t getQualityScoreSum() const
    {
        const Read &read = fragment_.getRead();
        return std::accumulate(read.getForwardQuality().begin(), read.getForwardQuality().end(), 0);
    }

    bool isAligned() const
    {
        return fragment_.isAligned();
    }

    bool isUniquelyAligned() const
    {
        return /*!fragment_.dodgy && */fragment_.isUniquelyAligned();
    }

    static uint64_t alignedBasesFromCigarOperation(const Cigar::value_type &cigarOperation)
    {
        const Cigar::Component decodedOperation = Cigar::decode(cigarOperation);
        return Cigar::ALIGN == decodedOperation.second ? decodedOperation.first : 0;
    }

    uint64_t getUniquelyAlignedBasesOutsideIndels() const
    {
        const alignment::Cigar::const_iterator cigarBegin = fragment_.cigarBuffer->begin() + fragment_.cigarOffset;
        const alignment::Cigar::const_iterator cigarEnd = cigarBegin + fragment_.cigarLength;

        return std::accumulate(cigarBegin, cigarEnd, 0,
                               bind(std::plus<uint64_t>(),
                                            _1,
                                            boost::bind(&FragmentMetadataTileStatsAdapter::alignedBasesFromCigarOperation, _2)));
    }

    uint64_t getMismatches() const
    {
        return fragment_.getMismatchCount();
    }

    uint64_t getEditDistance() const
    {
        return fragment_.getEditDistance();
    }

    unsigned getAlignmentScore() const
    {
        return fragment_.getAlignmentScore();
    }

    bool hasAlignmentScore() const
    {
        return fragment_.hasAlignmentScore();
    }

    std::vector<char>::const_iterator getForwardSequenceBegin() const
    {
        return fragment_.getRead().getForwardSequence().begin();
    }

    std::vector<char>::const_iterator getForwardSequenceEnd() const
    {
        return fragment_.getRead().getForwardSequence().end();
    }

    FragmentMetadata::ConstMismatchCycleIterator mismatchCyclesBegin() const
    {
        return fragment_.getMismatchCyclesBegin();
    }

    FragmentMetadata::ConstMismatchCycleIterator mismatchCyclesEnd() const
    {
        return fragment_.getMismatchCyclesEnd();
    }

    unsigned short getAdapterClipped() const
    {
        return fragment_.adapterClipped_;
    }

    friend std::ostream &operator << (std::ostream &os, const FragmentMetadataTileStatsAdapter &fragment)
    {
        return os << fragment.fragment_;
    }
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_FRAGMENT_METADATA_TILE_STATS_ADAPTER_HH
