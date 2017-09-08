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
 ** \file BinSorter.cpp
 **
 ** Reorders aligned data and stores results in bam file.
 ** 
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "bam/Bam.hh"

#include "build/BinSorter.hh"
#include "common/Memory.hh"

namespace isaac
{
namespace build
{

uint64_t BinSorter::serialize(
    BinData &binData,
    boost::ptr_vector<boost::iostreams::filtering_ostream> &bgzfStreams,
    boost::ptr_vector<bam::BamIndexPart> &bamIndexParts)
{
    if (!binData.getUniqueRecordsCount())
    {
        return 0;
    }
    ISAAC_THREAD_CERR << "Sorting offsets for bam " << binData.bin_ << std::endl;

    bamSerializer_.prepareForBam(contigLists_.front(), binData.data_, binData, binData.additionalCigars_, binData.splitInfoList_);

    ISAAC_THREAD_CERR << "Sorting offsets for bam done " << binData.bin_ << std::endl;

    ISAAC_THREAD_CERR << "Serializing records: " << binData.getUniqueRecordsCount() <<  " of them for bin " << binData.bin_ << std::endl;

    std::time_t serTimeStart = common::time();

    if (binData.isUnalignedBin())
    {
        uint64_t offset = 0;
        while(binData.data_.size() != offset)
        {
            const io::FragmentAccessor &fragment = binData.data_.getFragment(offset);
//            ISAAC_THREAD_CERR << "storeUnaligned: " << offset << "/" << binData.data_.size() << " " << fragment << std::endl;
            bamSerializer_.storeUnaligned(fragment, bgzfStreams, bamIndexParts, binData.bamAdapter_(fragment));
            offset += fragment.getTotalLength();
        }
    }
    else
    {
        BOOST_FOREACH(const PackedFragmentBuffer::Index& idx, binData)
        {
            // realigning reads that don't belong to the bin is not very useful
            // also, it can move the read position and cause more than one copy of the
            // read to be stored in the bam file.
            if (binData.bin_.hasPosition(idx.pos_))
            {
                const io::FragmentAccessor &fragment = binData.data_.getFragment(idx);
                bamSerializer_.storeAligned(fragment, bgzfStreams, bamIndexParts, binData.bamAdapter_(idx, fragment));
            }
            //else the fragment got split into a bit that does not belong to the current bin. it will get stored by another bin BinSorter.
        }
    }

    BOOST_FOREACH(boost::iostreams::filtering_ostream &bgzfStream, bgzfStreams)
    {
        ISAAC_VERIFY_MSG(bgzfStream.strict_sync(), "Expecting the compressor to flush all the data");
    }

    std::time_t serTimeEnd = common::time();

    ISAAC_THREAD_CERR << "Serializing records done: " << binData.getUniqueRecordsCount() <<  " of them for bin " << binData.bin_ << " in " << ::std::difftime(serTimeEnd, serTimeStart) << "seconds." << std::endl;
    return binData.size();
}

void BinSorter::resolveDuplicates(
    BinData &binData,
    BuildStats &buildStats)
{
    ISAAC_THREAD_CERR << "Resolving duplicates for bin " << binData.bin_ << std::endl;

    NotAFilter().filterInput(binData.data_, binData.seIdx_.begin(), binData.seIdx_.end(), buildStats, binData.binStatsIndex_, std::back_inserter(binData));
    if (keepDuplicates_ && !markDuplicates_)
    {
        NotAFilter().filterInput(binData.data_, binData.rIdx_.begin(), binData.rIdx_.end(), buildStats, binData.binStatsIndex_, std::back_inserter(binData));
        NotAFilter().filterInput(binData.data_, binData.fIdx_.begin(), binData.fIdx_.end(), buildStats, binData.binStatsIndex_, std::back_inserter(binData));
    }
    else
    {
        if (singleLibrarySamples_)
        {
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                RSDuplicateFilter<true>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.rIdx_.begin(), binData.rIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                FDuplicateFilter<true>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.fIdx_.begin(), binData.fIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
        }
        else
        {
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                RSDuplicateFilter<false>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.rIdx_.begin(), binData.rIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
            DuplicatePairEndFilter(keepDuplicates_).filterInput(
                FDuplicateFilter<false>(binData.barcodeBamMapping_.getSampleIndexMap()),
                binData.data_, binData.fIdx_.begin(), binData.fIdx_.end(),
                buildStats, binData.binStatsIndex_, std::back_inserter(binData));
        }
    }

    // we will not be needing these anymore. Free up some memory so that other bins get a chance to start earlier
    binData.unreserveIndexes();

    ISAAC_THREAD_CERR << "Resolving duplicates done for bin " << binData.bin_ << std::endl;
}

} // namespace build
} // namespace isaac
