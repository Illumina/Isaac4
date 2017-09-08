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
 ** \file GappedAligner.cpp
 **
 ** \brief See GappedAligner.hh
 ** 
 ** \author Come Raczy
 **/
#include "alignment/templateBuilder/GappedAligner.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

const unsigned short GappedAligner::UNINITIALIZED_OFFSET_MAGIC;
const unsigned short GappedAligner::REPEAT_OFFSET_MAGIC;

static std::vector<char> initVector;

GappedAligner::GappedAligner(
    const bool collectMismatchCycles,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const bool smartSmithWaterman,
    const unsigned smithWatermanGapSizeMax,
    const AlignmentCfg &alignmentCfg)
    : AlignerBase(collectMismatchCycles, alignmentCfg)
    , smartSmithWaterman_(smartSmithWaterman)
    , smithWatermanGapSizeMax_(smithWatermanGapSizeMax)
    , bandedSmithWaterman16_(alignmentCfg.gapMatchScore_, alignmentCfg.gapMismatchScore_, -alignmentCfg.gapOpenScore_, -alignmentCfg.gapExtendScore_,
                           flowcell::getMaxTotalReadLength(flowcellLayoutList))
    , bandedSmithWaterman32_(alignmentCfg.gapMatchScore_, alignmentCfg.gapMismatchScore_, -alignmentCfg.gapOpenScore_, -alignmentCfg.gapExtendScore_,
                       flowcell::getMaxTotalReadLength(flowcellLayoutList))
    , bandedSmithWaterman64_(alignmentCfg.gapMatchScore_, alignmentCfg.gapMismatchScore_, -alignmentCfg.gapOpenScore_, -alignmentCfg.gapExtendScore_,
                       flowcell::getMaxTotalReadLength(flowcellLayoutList))
    , hashedQueryTile_(2, -1U)
    , hashedQueryCluster_(2, -1U)
    , hashedQueryReadIndex_(2, -1U)
    , hashedQueryBegin_(2, initVector.begin())
    , hashedQueryEnd_(2, initVector.begin())
    , trackedOffsets_(QUERY_LENGTH_MAX, 0)
{
    queryKmerOffsets_[0].resize(oligo::MaxKmer<HASH_KMER_LENGTH, unsigned short>::value + 1, UNINITIALIZED_OFFSET_MAGIC);
    queryKmerOffsets_[1].resize(oligo::MaxKmer<HASH_KMER_LENGTH, unsigned short>::value + 1, UNINITIALIZED_OFFSET_MAGIC);
}

/// calculate the left and right flanks of the database WRT the query
static std::pair<unsigned, unsigned> getFlanks(
    const int64_t strandPosition,
    const unsigned readLength,
    const uint64_t referenceSize,
    const unsigned widestGapSize)
{
    assert(widestGapSize >= 2);
    if (strandPosition >= widestGapSize/2)
    {
        // take into account the rounding
        if (strandPosition + readLength + (widestGapSize - widestGapSize/2) < static_cast<int64_t>(referenceSize))
        {
            const unsigned left = widestGapSize/2;
            const unsigned right = widestGapSize - left - 1;
            return std::pair<unsigned, unsigned>(left, right);
        }
        else
        {
            assert((int64_t)referenceSize >= readLength + strandPosition);
            const unsigned right = referenceSize - readLength - strandPosition;
            const unsigned left = widestGapSize - right - 1;
            return std::pair<unsigned, unsigned>(left, right);
        }
    }
    else
    {
        assert(strandPosition >= 0);
        const unsigned left = strandPosition;
        const unsigned right = widestGapSize - left - 1;
        return std::pair<unsigned, unsigned>(left, right);
    }
}

/**
 * \brief as smith-waterman takes about 50% of the time with only 2% of results accepted, it is worth
 *        having a cheap way of checking whether the expensive gapped alignment is going to provide a meaningful outcome
 */
bool GappedAligner::makesSenseToGapAlign(
    const unsigned tile, const unsigned cluster, const unsigned read, const bool reverse,
    const std::vector<char>::const_iterator queryBegin,
    const std::vector<char>::const_iterator queryEnd,
    const reference::Contig::const_iterator databaseBegin,
    const reference::Contig::const_iterator databaseEnd)
{
    if (hashedQueryTile_[reverse] != tile || hashedQueryCluster_[reverse] != cluster ||
        hashedQueryReadIndex_[reverse] != read || hashedQueryBegin_[reverse] != queryBegin || hashedQueryEnd_[reverse] != queryEnd)
    {
        oligo::KmerGenerator<HASH_KMER_LENGTH, unsigned, std::vector<char>::const_iterator> queryKmerGenerator(queryBegin, queryEnd);
        std::fill(queryKmerOffsets_[reverse].begin(), queryKmerOffsets_[reverse].end(), UNINITIALIZED_OFFSET_MAGIC);
        unsigned kmer;
        std::vector<char>::const_iterator queryIt;
        while (queryKmerGenerator.next(kmer, queryIt))
        {
            if (UNINITIALIZED_OFFSET_MAGIC == queryKmerOffsets_[reverse][kmer])
            {
                queryKmerOffsets_[reverse][kmer] = (queryIt - queryBegin);
            }
            else
            {
                queryKmerOffsets_[reverse][kmer] = REPEAT_OFFSET_MAGIC;
            }
        }
        hashedQueryTile_[reverse] = tile;
        hashedQueryCluster_[reverse] = cluster;
        hashedQueryReadIndex_[reverse] = read;
        hashedQueryBegin_[reverse] = queryBegin;
        hashedQueryEnd_[reverse] = queryEnd;
    }

    const int queryLength = std::distance(queryBegin, queryEnd);
    ISAAC_ASSERT_MSG((queryLength * 2 + std::size_t(std::distance(databaseBegin, databaseEnd))) < trackedOffsets_.size(),
                     "number of tracked offsets is too small: " << trackedOffsets_.size() << " required: " <<
                     (queryLength * 2 + std::distance(databaseBegin, databaseEnd)));
    std::fill(trackedOffsets_.begin(), trackedOffsets_.begin() + queryLength * 2 + std::size_t(std::distance(databaseBegin, databaseEnd)), 0);

    oligo::KmerGenerator<HASH_KMER_LENGTH, unsigned, reference::Contig::const_iterator> databaseKmerGenerator(databaseBegin, databaseEnd);
    reference::Contig::const_iterator databaseIt;
    int lastConfirmedOffset = std::numeric_limits<int>::max();
    unsigned kmer;
    while (databaseKmerGenerator.next(kmer, databaseIt))
    {
        const int queryOffset = queryKmerOffsets_[reverse][kmer];
        if (REPEAT_OFFSET_MAGIC == queryOffset)
        {
            // Ambiguous mapping between reference and data,
            // Ignore, assume this is not informative...
        }
        else if (UNINITIALIZED_OFFSET_MAGIC != queryOffset)
        {
            const int databaseOffset = std::distance(databaseBegin, databaseIt);
            const int firstBaseOffset = databaseOffset - queryOffset + queryLength;
            // look for two sufficiently confirmed anchoring points that will disagree about the read offset.
            // If not found, sequence is unlikely to produce gaps with smith-waterman.
            if (++trackedOffsets_[firstBaseOffset] == SUFFICIENT_NUMBER_OF_HITS)
            {
                if (std::numeric_limits<int>::max() == lastConfirmedOffset)
                {
                    lastConfirmedOffset = firstBaseOffset;
                }
                else if (lastConfirmedOffset != firstBaseOffset)
                {
                    return true;
                }
            }
        }
    }

    return false;
}

template <typename BswT>
unsigned GappedAligner::alignGapped(
    BswT &bandedSmithWaterman,
    const bool smartSmithWaterman,
    const flowcell::ReadMetadata &readMetadata,
    const templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    const reference::ContigList &contigList,
    FragmentMetadata &fragmentMetadata,
    Cigar &cigarBuffer)
{
    const unsigned cigarOffset = cigarBuffer.size();
    fragmentMetadata.resetAlignment();
    fragmentMetadata.resetClipping();

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignGapped: after resetClipping: " << fragmentMetadata);

    const Read &read = fragmentMetadata.getRead();
    const bool reverse = fragmentMetadata.reverse;
    const std::vector<char> &sequence = read.getStrandSequence(reverse);
    const reference::Contig &contig = contigList[fragmentMetadata.contigId];

    std::vector<char>::const_iterator sequenceBegin = sequence.begin();
    std::vector<char>::const_iterator sequenceEnd = sequence.end();

    adapterClipper.clip(contig, fragmentMetadata, sequenceBegin, sequenceEnd);

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignGapped: after adapterClipper.clip: " << fragmentMetadata);
    clipReadMasking(read, fragmentMetadata, sequenceBegin, sequenceEnd);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignGapped: after clipReadMasking: " << fragmentMetadata);

    clipReference(contig.size(), fragmentMetadata.position, sequenceBegin, sequenceEnd);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignGapped: after clipReference: " << fragmentMetadata);

    const unsigned firstMappedBaseOffset = std::distance(sequence.begin(), sequenceBegin);
    if (firstMappedBaseOffset)
    {
        cigarBuffer.addOperation(firstMappedBaseOffset, Cigar::SOFT_CLIP);
    }

    const unsigned sequenceLength = std::distance(sequenceBegin, sequenceEnd);

    // position of the fragment on the strand
    int64_t strandPosition = fragmentMetadata.position;
    ISAAC_ASSERT_MSG(0 <= strandPosition, "alignUngapped should have clipped reads beginning before the reference");

    // no gapped alignment if the reference is too short
    if (static_cast<int64_t>(contig.size()) < sequenceLength + strandPosition + bandedSmithWaterman.WIDEST_GAP_SIZE)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "alignGapped: reference too short!");
        return 0;
    }
    // find appropriate beginning and end for the database
    const std::pair<unsigned, unsigned> flanks = getFlanks(strandPosition, sequenceLength, contig.size(), bandedSmithWaterman.WIDEST_GAP_SIZE);
    assert(flanks.first + flanks.second == bandedSmithWaterman.WIDEST_GAP_SIZE - 1);
    assert(flanks.first <= strandPosition);
    assert(strandPosition + sequenceLength + flanks.second <= (int64_t)contig.size());
    const reference::Contig::const_iterator databaseBegin = contig.begin() + strandPosition - flanks.first;
    const reference::Contig::const_iterator databaseEnd = databaseBegin + flanks.first + sequenceLength + flanks.second;

    if (smartSmithWaterman && !makesSenseToGapAlign(
        fragmentMetadata.getCluster().getTile(), fragmentMetadata.getCluster().getId(),
        fragmentMetadata.getReadIndex(), fragmentMetadata.isReverse(),
        sequenceBegin, sequenceEnd, databaseBegin, databaseEnd))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "Gap-aligning does not make sense" << common::makeFastIoString(sequenceBegin, sequenceEnd) <<
            " against " << common::makeFastIoString(databaseBegin, databaseEnd));
        return 0;
    }


    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "Gap-aligning " << common::makeFastIoString(sequenceBegin, sequenceEnd) <<
        " against " << common::makeFastIoString(databaseBegin, databaseEnd) << " strandPosition:"<<strandPosition);
    strandPosition += bandedSmithWaterman.align(sequenceBegin, sequenceEnd, databaseBegin, databaseEnd, cigarBuffer);

    if (firstMappedBaseOffset)
    {
        const Cigar::Component firstComponent = Cigar::decode(cigarBuffer.at(cigarOffset + 1));
        // avoid two softclips in a row
        if (Cigar::SOFT_CLIP == firstComponent.second)
        {
            cigarBuffer.erase(cigarBuffer.begin() + cigarOffset);
            cigarBuffer.updateOperation(cigarOffset, firstMappedBaseOffset + firstComponent.first, Cigar::SOFT_CLIP);
        }
    }

    const unsigned clipEndBases = std::distance(sequenceEnd, sequence.end());
    if (clipEndBases)
    {
        const Cigar::Component lastComponent = Cigar::decode(cigarBuffer.back());
        if (Cigar::SOFT_CLIP == lastComponent.second)
        {
            cigarBuffer.updateOperation(cigarBuffer.size() - 1, clipEndBases + lastComponent.first, Cigar::SOFT_CLIP);
        }
        else
        {
            cigarBuffer.addOperation(clipEndBases, Cigar::SOFT_CLIP);
        }
    }

    // adjust the start position of the fragment
    strandPosition -= flanks.first;

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "gapped CIGAR: " <<
//                                           alignment::Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()) << " strandPosition:"<<strandPosition);

    const unsigned matchCount = updateFragmentCigar(
        readMetadata, contigList, fragmentMetadata,
        fragmentMetadata.reverse, fragmentMetadata.contigId, strandPosition, cigarBuffer, cigarOffset);

    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "gapped: " << fragmentMetadata);

    return matchCount;
}

template <typename BswT>
bool GappedAligner::realignOne(
    BswT &bandedSmithWaterman,
    const bool smartSmithWaterman,
    const unsigned smitWatermanGapsMax,
    const reference::ContigList& contigList,
    const flowcell::ReadMetadata& readMetadata,
    FragmentMetadata& fragmentMetadata,
    templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
    Cigar& cigarBuffer)
{
    adapterClipper.checkInitStrand(fragmentMetadata, contigList[fragmentMetadata.contigId]);
    FragmentMetadata old = fragmentMetadata;
    const unsigned matchCount = alignGapped(smartSmithWaterman, readMetadata, adapterClipper, contigList, fragmentMetadata, cigarBuffer);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "    Gap-aligned: " << fragmentMetadata);
    if (matchCount &&
        // make sure we don't accept sw just moving ungapped alignments around. It only confuses the high-level logic
        fragmentMetadata.gapCount && fragmentMetadata.gapCount <= smitWatermanGapsMax &&
        fragmentMetadata.isBetterGapped(old))
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "    Using gap-aligned: " << fragmentMetadata);
        return true;
    }
    return false;
}

/**
 * \brief Will realign all bad ungapped alignments. If smart filtering is enabled, will realign first one regardless
 */
template <typename BswT>
bool GappedAligner::realignBadUngappedAlignments(
    BswT &bandedSmithWaterman,
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    FragmentMetadataList &fragments,
    templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
    Cigar &cigarBuffer)
{
    bool gappedFound = false;
    bool first = true;
    BOOST_FOREACH (const FragmentMetadata &fragmentMetadata, std::make_pair(fragments.begin(), fragments.end()))
    {
        // don't realign those that already have gaps detected by other means
        if (!fragmentMetadata.decoyAlignment && !fragmentMetadata.gapCount)
        {

            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(fragmentMetadata.getCluster().getId(), "    Original    : " << fragmentMetadata);
            if (bandedSmithWaterman.mismatchesMin_ <= fragmentMetadata.mismatchCount)
            {
                FragmentMetadata tmp = fragmentMetadata;
                if (realignOne(bandedSmithWaterman,
                    smartSmithWaterman_ && !first, smitWatermanGapsMax, contigList, readMetadata, tmp, adapterClipper, cigarBuffer))
                {
                    ISAAC_ASSERT_MSG(fragments.size() != fragments.capacity(), "Out of capacity in realignBadUngappedAlignments:" << fragments.capacity());
                    fragments.push_back(tmp);
                    gappedFound = true;
                }
            }
        }
        first = false;
    }

    if (gappedFound)
    {
        // gapped alignment and adapter trimming may have adjusted the alignment position
        std::sort(fragments.begin(), fragments.end());
        fragments.erase(std::unique(fragments.begin(), fragments.end()), fragments.end());
        putBestOnTop<true>(fragments);
        return true;
    }
    return false;
}

unsigned GappedAligner::alignGapped(
        const bool smartSmithWaterman,
        const flowcell::ReadMetadata &readMetadata,
        const FragmentSequencingAdapterClipper &adapterClipper,
        const reference::ContigList &contigList,
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer)
{
    switch(smithWatermanGapSizeMax_)
    {
    case 16:
        return alignGapped(
            bandedSmithWaterman16_,
            smartSmithWaterman, readMetadata, adapterClipper, contigList, fragmentMetadata, cigarBuffer);
    case 32:
        return alignGapped(
            bandedSmithWaterman32_,
            smartSmithWaterman, readMetadata, adapterClipper, contigList, fragmentMetadata, cigarBuffer);
    case 64:
        return alignGapped(
            bandedSmithWaterman64_,
            smartSmithWaterman, readMetadata, adapterClipper, contigList, fragmentMetadata, cigarBuffer);
    default:
        BOOST_THROW_EXCEPTION(common::InvalidParameterException("Unsupported smithWatermanGapSizeMax"));
    }

    // should not get here
    ISAAC_ASSERT_MSG(false, "Can't be here")
    return 0;
}

bool GappedAligner::realignBadUngappedAlignments(
    const unsigned gappedMismatchesMax,
    const unsigned smitWatermanGapsMax,
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata &readMetadata,
    FragmentMetadataList &fragmentList,
    FragmentSequencingAdapterClipper &adapterClipper,
    Cigar &cigarBuffer)
{
    switch(smithWatermanGapSizeMax_)
    {
    case 16:
        return realignBadUngappedAlignments(
            bandedSmithWaterman16_,
            gappedMismatchesMax, smitWatermanGapsMax, contigList, readMetadata, fragmentList, adapterClipper, cigarBuffer);
    case 32:
        return realignBadUngappedAlignments(
            bandedSmithWaterman32_,
            gappedMismatchesMax, smitWatermanGapsMax, contigList, readMetadata, fragmentList, adapterClipper, cigarBuffer);
    case 64:
        return realignBadUngappedAlignments(
            bandedSmithWaterman64_,
            gappedMismatchesMax, smitWatermanGapsMax, contigList, readMetadata, fragmentList, adapterClipper, cigarBuffer);
    default:
        BOOST_THROW_EXCEPTION(common::InvalidParameterException("Unsupported smithWatermanGapSizeMax"));
    }
    // should not get here
    ISAAC_ASSERT_MSG(false, "Can't be here")
    return false;
}


} // namespace templateBuilder
} // namespace alignment
} // namespace isaac
