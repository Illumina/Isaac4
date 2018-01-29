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
 ** \file FindHashMatchesTransition.cpp
 **
 ** \brief see FindHashMatchesTransition.hh
 **
 ** \author Roman Petrovski
 **/

#include "flowcell/Layout.hh"
#include "alignment/HashMatchFinder.hh"
#include "alignment/Quality.hh"
#include "oligo/KmerGenerator.hpp"
#include "reference/Seed.hh"

namespace isaac
{
namespace alignment
{

template <typename ReferenceHash>
SeedHashMatchFinder<ReferenceHash>::SeedHashMatchFinder(
    const ReferenceHash& referenceHash,
    const unsigned seedBaseQualityMin)
    : referenceHash_(referenceHash)
    , seedBaseQualityMin_(seedBaseQualityMin)
{
}

inline reference::ContigList::Offset getForwardAlignmentOffset(
    const int seedOffset,
    const reference::ContigList::Offset seedLocation)
{
    return seedLocation - seedOffset;
}

template <typename SeedT>
inline reference::ContigList::Offset getReverseAlignmentOffset(
    const unsigned seedOffset, const unsigned readLength,
    const reference::ContigList::Offset seedLocation)
{
    const unsigned seedLength = SeedT::SEED_LENGTH;
    const unsigned seedStep = SeedT::STEP;
    return seedLocation + seedLength + seedOffset - readLength +
        //step-1 seed reverse match is 0 bases off
        //step-2 seed reverse match is 1 base off
        //step-3 seed reverse match is 2 bases off
        //...
        (1 - seedStep);
}

template <typename ReferenceHash>
template <bool reverse, bool filterContigs>
std::size_t iSAAC_PROFILING_NOINLINE SeedHashMatchFinder<ReferenceHash>::getSeedAlignmentPositions(
    const reference::ContigList &contigList,
    const unsigned seedOffset,
    const flowcell::ReadMetadata &readMetadata,
    const unsigned filterContigId,
    typename ReferenceHash::MatchRange matchPositions,
    ReferenceOffsetList& referenceOffsets) const
{
    const std::size_t before = referenceOffsets.size();

//    ISAAC_ASSERT_MSG(matchPositions.second == std::adjacent_find(matchPositions.first, matchPositions.second),
//                     "Duplicate matches unexpected:" << *std::adjacent_find(matchPositions.first, matchPositions.second));
    if (filterContigs)
    {
        const reference::ContigList::Offset filterOffsetBegin = contigList.beginOffset(filterContigId);
        const reference::ContigList::Offset filterOffsetEnd = contigList.endOffset(filterContigId);
//        ISAAC_THREAD_CERR << " filterOffsetBegin:" << filterOffsetBegin << " filterOffsetBegin:" << filterOffsetBegin << std::endl;
        for(auto it = matchPositions.first; it != matchPositions.second; ++it)
        {
            const reference::ContigList::Offset &referenceOffset = *it;
            if (filterOffsetBegin <= referenceOffset && referenceOffset < filterOffsetEnd)
            {
                if (reverse)
                {
                    referenceOffsets.push_back(getReverseAlignmentOffset<reference::Seed<KmerT> >(seedOffset, readMetadata.getLength(), referenceOffset));
//                    ISAAC_THREAD_CERR << "rv: " << " seedOffset:" << seedOffset << " referenceOffset:" << referenceOffset << " pos:" << referenceOffsets.back() << std::endl;
                }
                else
                {
                    referenceOffsets.push_back(getForwardAlignmentOffset(seedOffset, referenceOffset));
//                    ISAAC_THREAD_CERR << "fw: " << " seedOffset:" << seedOffset << " referenceOffset:" << referenceOffset << " pos:" << referenceOffsets.back() << std::endl;
                }
            }
        }
    }
    else
    {
        for(auto it = matchPositions.first; it != matchPositions.second; ++it)
        {
            const reference::ContigList::Offset &referenceOffset = *it;
            if (reverse)
            {
                referenceOffsets.push_back(getReverseAlignmentOffset<reference::Seed<KmerT> >(seedOffset, readMetadata.getLength(), referenceOffset));
//                ISAAC_THREAD_CERR << "rv: " << " seedOffset:" << seedOffset << " referenceOffset:" << referenceOffset << " pos:" << referenceOffsets.back() << std::endl;
            }
            else
            {
                referenceOffsets.push_back(getForwardAlignmentOffset(seedOffset, referenceOffset));
//                ISAAC_THREAD_CERR << "fw: " << " seedOffset:" << seedOffset << " referenceOffset:" << referenceOffset << " pos:" << referenceOffsets.back() << std::endl;
            }
        }
    }

    return referenceOffsets.size() - before;
}

template <typename ReferenceHash, unsigned seedsPerMatchMax>
ClusterHashMatchFinder<ReferenceHash, seedsPerMatchMax>::ClusterHashMatchFinder(
    const ReferenceHash& referenceHash,
    const std::size_t candidateMatchesMax,
    const unsigned seedBaseQualityMin,
    const unsigned seedRepeatsMax) :
    SeedHashMatchFinder<ReferenceHash>(referenceHash, seedBaseQualityMin),
    candidateMatchesMax_(candidateMatchesMax)
{
}


/**
 * \brief advances begin until begin == end || *begin == hit || *begin > hit. Inserts all *begin < hit into keep.
 *        Inserts *begin == hit into promoteBuffer
 *
 * \return true if match found, false otherwise.
 */
inline bool iSAAC_PROFILING_NOINLINE mergeSeedHit(
    ReferenceOffsetList::iterator &begin,
    const ReferenceOffsetList::iterator end,
    const reference::ContigList::Offset hit,
    ReferenceOffsetList::iterator &keep,
    ReferenceOffsetList &promoteBuffer)
{
    ISAAC_ASSERT_MSG(begin <= end, "std::distance(begin, end):" << std::distance(begin, end));
    for (;end != begin; ++begin)
    {
        if (*begin < hit)
        {
            *keep = *begin;
            ++keep;
        }
        else if (*begin == hit)
        {
            promoteBuffer.push_back(hit);
            ++begin;
            return true;
        }
        else
        {
            return false;
        }
    }
    return false;
}

struct Unmerged
{
    // cursor
    ReferenceOffsetList::iterator current_;
    // end of input
    const ReferenceOffsetList::iterator end_;
    // end of input once all the promoted have been removed and
    // all the ones that stayed got moved up to fill the hole.
    ReferenceOffsetList::iterator newEnd_;
};


template <unsigned maxSupportingSeeds>
void mergeSeedHiT(
    ReferenceOffset hit,
    Unmerged unmerged[],
    ReferenceOffsetLists& mergeBuffers)
{
    if (!mergeSeedHit(unmerged[maxSupportingSeeds - 1].current_, unmerged[maxSupportingSeeds - 1].end_, hit, unmerged[maxSupportingSeeds - 1].newEnd_, mergeBuffers[maxSupportingSeeds]))
    {
        mergeSeedHiT<maxSupportingSeeds - 1>(hit, unmerged, mergeBuffers);
    }
}

template <>
void mergeSeedHiT<1>(
    ReferenceOffset hit,
    Unmerged unmerged[],
    ReferenceOffsetLists& mergeBuffers)
{
    mergeBuffers[1].push_back(hit);
}

/**
 * \brief merge two sorted sequences into one. The input sequences are in the same buffer
 */
bool iSAAC_PROFILING_NOINLINE mergeUnmerged(
    Unmerged &unmerged, ReferenceOffsetList& mergeBuffer, ReferenceOffsetList& tmp)
{
    unmerged.newEnd_ = unmerged.current_ != unmerged.newEnd_ ? std::move(unmerged.current_, unmerged.end_, unmerged.newEnd_) : unmerged.end_;
    //        mergeBuffers[i].erase(, mergeBuffers[i+1].end());
    tmp.clear();
    std::merge(mergeBuffer.begin(), unmerged.newEnd_, unmerged.end_, mergeBuffer.end(), std::back_inserter(tmp));
    mergeBuffer.clear();
    mergeBuffer.insert(mergeBuffer.end(), tmp.begin(), tmp.end());
    tmp.clear();
    return !mergeBuffer.empty();
}

/**
 * \return the highest count of seeds supporting a match
 */
template <unsigned maxSupportingSeeds> unsigned mergeUnmerged(
    Unmerged *unmerged,
    ReferenceOffsetLists& mergeBuffers,
    ReferenceOffsetList& tmp)
{
    unsigned ret = mergeUnmerged(unmerged[maxSupportingSeeds], mergeBuffers[maxSupportingSeeds], tmp) ? maxSupportingSeeds : 0;
    unsigned next = mergeUnmerged<maxSupportingSeeds - 1>(unmerged, mergeBuffers, tmp);
    return ret ? ret : next;
}

template <> unsigned mergeUnmerged<0>(
    Unmerged *unmerged,
    ReferenceOffsetLists& mergeBuffers,
    ReferenceOffsetList& tmp)
{
    return 0;
}

// please implement a specialization for your maxSupportingSeeds
template <unsigned maxSupportingSeeds>
unsigned iSAAC_PROFILING_NOINLINE mergeSeedHits(ReferenceOffsetLists& mergeBuffers)
{
    ISAAC_ASSERT_MSG(false, "Implementation must be explicity specialized");
    return -1U;
}

template <>
unsigned iSAAC_PROFILING_NOINLINE mergeSeedHits<4>(ReferenceOffsetLists& mergeBuffers)
{
    static const unsigned maxSupportingSeeds = 4;
    const ReferenceOffsetList& seedHits = mergeBuffers[0];
    if (seedHits.empty())
    {
        return 0;
    }

    Unmerged unmerged[maxSupportingSeeds + 1] =
    {
        {mergeBuffers[1].begin(), mergeBuffers[1].end(), mergeBuffers[1].begin()}, //[0] is unused. Just a filler to avoid warnings on gcc 4.7
        {mergeBuffers[1].begin(), mergeBuffers[1].end(), mergeBuffers[1].begin()},
        {mergeBuffers[2].begin(), mergeBuffers[2].end(), mergeBuffers[2].begin()},
        {mergeBuffers[3].begin(), mergeBuffers[3].end(), mergeBuffers[3].begin()},
        {mergeBuffers[4].begin(), mergeBuffers[4].end(), mergeBuffers[4].begin()},
//        {mergeBuffers[5].begin(), mergeBuffers[5].end(), mergeBuffers[5].begin()},
//        {mergeBuffers[6].begin(), mergeBuffers[6].end(), mergeBuffers[6].begin()},
//        {mergeBuffers[7].begin(), mergeBuffers[7].end(), mergeBuffers[7].begin()},
//        {mergeBuffers[8].begin(), mergeBuffers[8].end(), mergeBuffers[8].begin()},
//        {mergeBuffers[9].begin(), mergeBuffers[9].end(), mergeBuffers[9].begin()},
//        {mergeBuffers[10].begin(), mergeBuffers[10].end(), mergeBuffers[10].begin()},
    };

    for (ReferenceOffset hit : seedHits)
    {
        mergeSeedHiT<maxSupportingSeeds>(hit, unmerged, mergeBuffers);
    }

    // done using input, clear buffer and pass it as tmp for mergeUnmerged
    mergeBuffers[0].clear();
    //WARNING! after mergeUnmerged the unmerged is not pointing correctly.
    return mergeUnmerged<maxSupportingSeeds>(unmerged, mergeBuffers, mergeBuffers[0]);

//    const bool have4 = mergeUnmerged(unmerged[4], mergeBuffers[4], mergeBuffers[0]);
//    const bool have3 = mergeUnmerged(unmerged[3], mergeBuffers[3], mergeBuffers[0]);
//
//    bool have2 = false;
//    bool have1 = false;
//    // performance HACK: we will stop when we have 4, so don't bother merging 2 and 1
//    if (!have4)
//    {
//        have2 = mergeUnmerged(unmerged[2], mergeBuffers[2], mergeBuffers[0]);
//        have1 = mergeUnmerged(unmerged[1], mergeBuffers[1], mergeBuffers[0]);
//    }
//
//    return have4 ? 4 : have3 ? 3 : have2 ? 2 : have1 ? 1 : 0;
}


template<bool reverse>
iSAAC_PROFILING_NOINLINE
void offsetsToMatches(const ReferenceOffsetList& mergeBuffer, Matches& matches)
{
    for (const ReferenceOffset offset : mergeBuffer)
    {
        matches.push_back(Match(offset, reverse));
    }
}

template <typename ReferenceHash, unsigned seedsPerMatchMax>
template<bool filterContigs, bool detectStructuralVariant>
void iSAAC_PROFILING_NOINLINE ClusterHashMatchFinder<ReferenceHash, seedsPerMatchMax>::buildMatchesIteratively(
    const reference::ContigList &contigList,
    const flowcell::ReadMetadata& readMetadata,
    const Cluster& cluster,
    const unsigned filterContigId,
    const SeedsHits& seedsHits,
    MatchLists& matchLists,
    ReferenceOffsetLists& fwMergeBuffers,
    ReferenceOffsetLists& rvMergeBuffers) const
{
    for (auto &buffer : fwMergeBuffers) {buffer.clear();}
    for (auto &buffer : rvMergeBuffers) {buffer.clear();}

    unsigned topSeeds = 0;
    for (const SeedHits &seedHits : seedsHits)
    {

        this->template getSeedAlignmentPositions<false, filterContigs>(
            contigList, seedHits.seedOffset_, readMetadata, filterContigId, seedHits.forwardMatches_, fwMergeBuffers[0]);
        topSeeds = std::max(topSeeds, mergeSeedHits<seedsPerMatchMax>(fwMergeBuffers));
        this->template getSeedAlignmentPositions<true, filterContigs>(
            contigList, seedHits.seedOffset_, readMetadata, filterContigId, seedHits.reverseMatches_, rvMergeBuffers[0]);
        topSeeds = std::max(topSeeds, mergeSeedHits<seedsPerMatchMax>(rvMergeBuffers));

        if (!detectStructuralVariant && LONG_READ_SEEDS_MIN <= topSeeds &&
            (fwMergeBuffers[topSeeds - 1].size() + rvMergeBuffers[topSeeds - 1].size()) &&
            (fwMergeBuffers[topSeeds].size() + rvMergeBuffers[topSeeds].size() +
            fwMergeBuffers[topSeeds - 1].size() + rvMergeBuffers[topSeeds - 1].size()) < candidateMatchesMax_)
        {
            break;
        }

        if (topSeeds == seedsPerMatchMax)
        {
            break;
        }
    }

    ISAAC_ASSERT_MSG(matchLists.size() > topSeeds, "Incorrectly allocated buffer. Expected: " << (topSeeds + 1) << "Got:" << matchLists.size());
    ISAAC_ASSERT_MSG(!seedsHits.empty(), "Unexpected empty seedsHits. Caller should ensure the input makes sense.");
    const unsigned minSeeds = std::min<unsigned>(seedsHits.size(), detectStructuralVariant ? SV_READ_SEEDS_MIN : SHORT_READ_SEEDS_MIN);

    if (minSeeds <= topSeeds)
    {
        const unsigned bestCandidateMatches = fwMergeBuffers[topSeeds].size() + rvMergeBuffers[topSeeds].size();
        const unsigned secondBestCandidateMatches = fwMergeBuffers[topSeeds - 1].size() + rvMergeBuffers[topSeeds - 1].size();
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "buildMatchesIteratively bestCandidateMatches: " << bestCandidateMatches);
        if ((bestCandidateMatches + secondBestCandidateMatches) < candidateMatchesMax_)
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "buildMatchesIteratively secondBestCandidateMatches:" << secondBestCandidateMatches);
            offsetsToMatches<false>(fwMergeBuffers[topSeeds], matchLists[topSeeds]);
            offsetsToMatches<true>(rvMergeBuffers[topSeeds], matchLists[topSeeds]);
//              if (!detectStructuralVariant)
            {
                offsetsToMatches<false>(fwMergeBuffers[topSeeds - 1], matchLists[topSeeds - 1]);
                offsetsToMatches<true>(rvMergeBuffers[topSeeds - 1], matchLists[topSeeds - 1]);
            }

            if (2 < topSeeds && !secondBestCandidateMatches)
            {
                const unsigned thirdBestCandidateMatches = fwMergeBuffers[topSeeds - 2].size() + rvMergeBuffers[topSeeds - 2].size();
                if (bestCandidateMatches + thirdBestCandidateMatches < candidateMatchesMax_)
                {
                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "buildMatchesIteratively thirdBestCandidateMatches:" << thirdBestCandidateMatches);

                    offsetsToMatches<false>(fwMergeBuffers[topSeeds - 2], matchLists[topSeeds - 2]);
                    offsetsToMatches<true>(rvMergeBuffers[topSeeds - 2], matchLists[topSeeds - 2]);

                    //This gives marginal improvement on scoring quality in simulations without any extra coverage.
                    //In fact the mapq 20+ coverage is reduced by a tiny bit due to higher accruacy scoring
//                    if (3 < topSeeds && !thirdBestCandidateMatches)
//                    {
//                        const unsigned fourthBestCandidateMatches = fwMergeBuffers[topSeeds - 3].size() + rvMergeBuffers[topSeeds - 3].size();
//                        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "buildMatchesIteratively fourthBestCandidateMatches:" << fourthBestCandidateMatches);
//                        if (bestCandidateMatches + fourthBestCandidateMatches < candidateMatchesMax_)
//                        {
//                            offsetsToMatches<false>(fwMergeBuffers[topSeeds - 3], matchLists[topSeeds - 3]);
//                            offsetsToMatches<true>(rvMergeBuffers[topSeeds - 3], matchLists[topSeeds - 3]);
//                        }
//                    }
                }
            }
        }
    }
    // else we either have too many candidate alignments or we have not used enough seeds to trust them.
}

template <typename ReferenceHash, unsigned seedsPerMatchMax>
std::size_t iSAAC_PROFILING_NOINLINE ClusterHashMatchFinder<ReferenceHash, seedsPerMatchMax>::collectSeedHits(
    const Cluster& cluster,
    const unsigned readIndex,
    const std::size_t seedRepeatThreshold,
    const unsigned endSeedOffset,
    SeedsHits& seedsHits) const
{
    struct BadBaseMasker
    {
        BadBaseMasker(const unsigned char seedBaseQualityMin = 0) : seedBaseQualityMin_(seedBaseQualityMin){}
        unsigned char seedBaseQualityMin_;
        unsigned char operator[](const char &base) const
        {
            return oligo::getQuality(base) < seedBaseQualityMin_ ?
                oligo::INVALID_OLIGO : static_cast<unsigned char>(base & oligo::BCL_BASE_MASK);
        }
    } translator(BaseT::seedBaseQualityMin_);

    typedef reference::Seed <KmerT> Seed;
    const BclClusters::const_iterator bclBegin = cluster.getBclData(readIndex);
    oligo::InterleavedKmerGenerator<Seed::KMER_BASES, typename Seed::KmerType, BclClusters::const_iterator, Seed::STEP, decltype(translator)>
        kmerGenerator(bclBegin, bclBegin + endSeedOffset, translator);

    std::size_t repeatSeeds = 0;
    KmerT seedKmer(0);
    BclClusters::const_iterator bclCurrent;
    while (kmerGenerator.next(seedKmer, bclCurrent))
    {
        const unsigned seedOffset = std::distance(bclBegin, bclCurrent);

        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(
            cluster.getId(), "seed at offset : " << seedOffset << " " <<
            (oligo::Bases<oligo::BITS_PER_BASE, KmerT>(seedKmer, oligo::KmerTraits<KmerT>::KMER_BASES)) << "/" <<
            (oligo::ReverseBases<oligo::BITS_PER_BASE, KmerT>(seedKmer, oligo::KmerTraits<KmerT>::KMER_BASES)) << " endSeedOffset:" << endSeedOffset);
        const typename ReferenceHash::MatchRange fwMatchRange = BaseT::referenceHash_.findMatches(seedKmer);
//        ISAAC_ASSERT_MSG(fwMatchRange.second == std::adjacent_find(fwMatchRange.first, fwMatchRange.second),
//                         "Duplicate matches unexpected:" << *std::adjacent_find(fwMatchRange.first, fwMatchRange.second) << " " << oligo::bases<2>(seedKmer, Seed::KMER_BASES));
//            for(auto it = fwMatchRange.first; it != fwMatchRange.second; ++it)
//            {
//                const reference::ContigList::Offset &referenceOffset = *it;
//                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "fw Hit offset:" << referenceOffset);
//            }
        if (std::size_t(std::distance(fwMatchRange.first, fwMatchRange.second)) >= seedRepeatThreshold)
        {
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "findReadMatches: " << seedOffset << " fwMatchRange: MatchRange(" << std::distance(fwMatchRange.first, fwMatchRange.second) << ")");
            ++repeatSeeds;
        }
        else
        {
            seedKmer = oligo::reverseComplement(seedKmer);

            const typename ReferenceHash::MatchRange rvMatchRange = BaseT::referenceHash_.findMatches(seedKmer);
//            ISAAC_ASSERT_MSG(rvMatchRange.second == std::adjacent_find(rvMatchRange.first, rvMatchRange.second),
//                             "Duplicate matches unexpected:" << *std::adjacent_find(rvMatchRange.first, rvMatchRange.second) << " " << oligo::bases<2>(seedKmer, Seed::KMER_BASES));
//            for(auto it = rvMatchRange.first; it != rvMatchRange.second; ++it)
//            {
//                const reference::ContigList::Offset &referenceOffset = *it;
//                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "rv Hit offset:" << referenceOffset);
//            }
            const SeedHits hits = { seedOffset, fwMatchRange, rvMatchRange };
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(cluster.getId(), "findReadMatches: " << seedOffset << " " << hits);
            if (hits.hitCount() >= seedRepeatThreshold)
            {
                ++repeatSeeds;
            }
            else if (!hits.empty()) //empty hits are either due to no match (unlikely in human) or seed base quality filtering.
            {
                seedsHits.push_back(hits);
                kmerGenerator.skip(KmerT::KMER_BASES - 1);
            }
        }
    }
    return repeatSeeds;
}

template <typename ReferenceHash, unsigned seedsPerMatchMax>
iSAAC_PROFILING_NOINLINE
std::size_t ClusterHashMatchFinder<ReferenceHash, seedsPerMatchMax>::findReadMatches(
    const reference::ContigList &contigList,
    const Cluster& cluster,
    const flowcell::ReadMetadata& readMetadata,
    const std::size_t seedRepeatThreshold,
    MatchLists& matchLists,
    ReferenceOffsetLists& fwMergeBuffers,
    ReferenceOffsetLists& rvMergeBuffers) const
{
    //ISAAC_ASSERT_CERR << "findReadMatches" << std::endl;
    ISAAC_ASSERT_MSG(matchLists.capacity() > seedsPerMatchMax,
                     "Insufficient capacity in matchLists:" << matchLists.capacity() << " for:" << seedsPerMatchMax << " seedsPerMatchMax");
    for (Matches &matches : matchLists) {matches.clear();}

    SeedsHits seedsHits;
    const std::size_t repeatSeeds = collectSeedHits(cluster, readMetadata.getIndex(), seedRepeatThreshold, readMetadata.getLength(), seedsHits);

    // demand LONG_READ_SEEDS_MIN unless read is too short, otherwise demand SHORT_READ_SEEDS_MIN.
    const unsigned seedsMin = std::min(LONG_READ_SEEDS_MIN, std::max(SHORT_READ_SEEDS_MIN, readMetadata.getLength() / 2 / SEED_LENGTH));
    if (seedsMin <= seedsHits.size())
    {
        std::sort(seedsHits.begin(), seedsHits.end());
        buildMatchesIteratively<false, false>(
            contigList, readMetadata, cluster, BaseT::NO_CONTIG_FILTER, seedsHits, matchLists, fwMergeBuffers, rvMergeBuffers);
    }

    return repeatSeeds;

}

template <typename ReferenceHash, unsigned seedsPerMatchMax>
std::size_t ClusterHashMatchFinder<ReferenceHash, seedsPerMatchMax>::findHeadAnchoredReadMatches(
    const reference::ContigList &contigList,
    const Cluster& cluster,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned filterContigId,
    const std::size_t seedRepeatThreshold,
    const unsigned headSeedOffsetMax,
    MatchLists& matchLists,
    ReferenceOffsetLists& fwMergeBuffers,
    ReferenceOffsetLists& rvMergeBuffers) const
{
    //ISAAC_ASSERT_CERR << "findReadMatches" << std::endl;
    ISAAC_ASSERT_MSG(matchLists.capacity() > seedsPerMatchMax,
                     "Insufficient capacity in matchLists:" << matchLists.capacity() << " for:" << seedsPerMatchMax << " seedsPerMatchMax");
    for (Matches &matches : matchLists) {matches.clear();}

    ISAAC_ASSERT_MSG(headSeedOffsetMax + reference::Seed<KmerT>::SEED_LENGTH <= readMetadata.getLength(), "Invalid headSeedOffsetMax:" << headSeedOffsetMax << " " << cluster);
//    if (headSeedOffsetMax + reference::Seed<KmerT>::SEED_LENGTH <= readMetadata.getLength())
    {
        SeedsHits seedsHits;
        const std::size_t repeatSeeds = collectSeedHits(cluster, readMetadata.getIndex(), seedRepeatThreshold, headSeedOffsetMax, seedsHits);
        if (SV_READ_SEEDS_MIN <= seedsHits.size())
        {
            std::sort(seedsHits.begin(), seedsHits.end());

            if (BaseT::NO_CONTIG_FILTER == filterContigId)
            {
                buildMatchesIteratively<false, true>(
                    contigList, readMetadata, cluster, BaseT::NO_CONTIG_FILTER, seedsHits, matchLists, fwMergeBuffers, rvMergeBuffers);
            }
            else
            {
                buildMatchesIteratively<true, true>(
                    contigList, readMetadata, cluster, filterContigId, seedsHits, matchLists, fwMergeBuffers, rvMergeBuffers);
            }
        }

        return repeatSeeds;
    }
    return 0;
}

template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::VeryShortKmerType>, 4>;

template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<10>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<11>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<12>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<13>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<14>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<15>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<16>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<17>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<18>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<19>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ClusterHashMatchFinder<reference::ReferenceHash<oligo::BasicKmerType<20>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;

} // namespace alignment
} // namespace isaac
