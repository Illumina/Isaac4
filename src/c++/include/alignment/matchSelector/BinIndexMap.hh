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
 ** \file BinIndexMap.hh
 **
 ** \brief Fragment buffer flushing and output file management.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH

#include <boost/foreach.hpp>

#include "alignment/MatchDistribution.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

class BinIndexMap: public std::vector<std::vector<unsigned> >
{
    /// the binSize from the MatchDistribution
    const unsigned binLength_;
public:
    BinIndexMap(
        const isaac::reference::SortedReferenceMetadata &sortedReference,
        const unsigned binLength)
        : binLength_(binLength)
    {
        // first push the contig and bin for unaligned clusters
        push_back(std::vector<unsigned>(1, 0));

        size_t currentBinIndex = 1;
        // now put in all the contig bins
        for (const isaac::reference::SortedReferenceMetadata::Contig &contig : sortedReference.getContigs())
        {
            push_back(std::vector<unsigned>((contig.totalBases_ + binLength_ - 1) / binLength_, 0));
            std::iota(back().begin(), back().end(), currentBinIndex);
            currentBinIndex = back().back() + 1;
        }
    }

    /**
     ** \brief convert a reference position on a contig into a bin index that
     ** can be used to identify either the file path or the stream associated
     ** to the ReferencePosition.
     **/
    std::size_t getBinIndex(const isaac::reference::ReferencePosition &referencePosition) const
    {
        const uint64_t contigId = referencePosition.getContigId();
        const std::vector<unsigned> &binIndexList = at(contigId + 1);
        const uint64_t position = referencePosition.getPosition();
        const uint64_t index = position / binLength_;
        ISAAC_ASSERT_MSG(binIndexList.size() > index, index << " required for vector of size " << back().size() << " " << referencePosition);
        return binIndexList[index];
    }

    /**
     * \return The first reference position that can be found in the bin
     */
    isaac::reference::ReferencePosition getBinFirstPos(const unsigned bin) const
    {
        std::vector<unsigned>::const_reference (std::vector<unsigned>::*f)() const = &std::vector<unsigned>::front;

        // use upper_bound to skip all the contigs that were so empty that they did not get mapped to a bin
        const_iterator binContigIterator = std::upper_bound(begin(), end(), bin, boost::bind(f, _2) > _1);
        ISAAC_ASSERT_MSG(begin() != binContigIterator, "Bin number has to be one of those we have a contig for: " << bin);
        //take a step back as we just skipped the last one we were looking for
        --binContigIterator;

        std::vector<unsigned>::const_iterator contigBinIterator =
            lower_bound(binContigIterator->begin(), binContigIterator->end(), bin);

        ISAAC_ASSERT_MSG(binContigIterator->end() != contigBinIterator, "Bin number must be present in the contig bins");

        const uint64_t contigId = binContigIterator - begin() - 1;
        const uint64_t position = binLength_ * (contigBinIterator - binContigIterator->begin());

        return isaac::reference::ReferencePosition(contigId, position);
    }

    /**
     * \return The first reference position that belongs to the subsequent bin. NOTE: for last bin in the contig
     *         there is not guarantee that no alignments will exist at this position and beyond. However, the amount
     *         of data aligning there should be considered minor and belonging to the last bin.
     */
    isaac::reference::ReferencePosition getBinFirstInvalidPos(const unsigned bin) const
    {
        std::vector<unsigned>::const_reference (std::vector<unsigned>::*f)() const = &std::vector<unsigned>::front;

        // use upper_bound to skip all the contigs that were so empty that they did not get mapped to a bin
        const_iterator binContigIterator = std::upper_bound(begin(), end(), bin, boost::bind(f, _2) > _1);
        ISAAC_ASSERT_MSG(begin() != binContigIterator, "Bin number has to be one of those we have a contig for: " << bin);
        //take a step back as we just skipped the last one we were looking for
        --binContigIterator;

        std::vector<unsigned>::const_iterator contigNextBinIterator =
            upper_bound(binContigIterator->begin(), binContigIterator->end(), bin);

        const uint64_t contigId = binContigIterator - begin() - 1;
        const uint64_t position = binLength_ * (contigNextBinIterator - binContigIterator->begin());

        return isaac::reference::ReferencePosition(contigId, position);
    }

    /**
     * \return The highest bin index to which the mapping is stored. Notice that there might be no bin with
     *         this index as it could have had no matches.
     */
    unsigned getHighestBinIndex() const
    {
        return back().back();
    }

    unsigned getTotalBins() const
    {
        return getHighestBinIndex() + 1;
    }

    unsigned getBinLength() const
    {
        return binLength_;
    }

    friend std::ostream& operator << (std::ostream& os, const BinIndexMap &binIndexMap)
    {
        return os << "BinIndexMap(" << binIndexMap.binLength_ <<"bl)";
    }
};

//
//inline std::ostream& operator << (std::ostream& os, const BinIndexMap &binIndexMap)
//{
//    std::string message = "\n";
//    unsigned index = 0;
//    BOOST_FOREACH(const std::vector<unsigned> &contigIndexList, binIndexMap)
//    {
//        message += boost::lexical_cast<std::string>(index) + ":";
//        message += contigIndexList.empty() ?
//            std::string("Empty bin list") :
//            (boost::format("%d bin indexes from %d to %d:") % contigIndexList.size() % contigIndexList.front() % contigIndexList.back()).str();
///*
//        BOOST_FOREACH(unsigned index, contigIndexList)
//        {
//            message += (boost::format(" %d") % index).str();
//        }
//*/
//        message += "\n";
//        ++index;
//    }
//    os << message;
//    return os;
//}

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_BIN_INDEX_MAP_HH
