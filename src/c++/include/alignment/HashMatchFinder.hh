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
 ** \file FindHashMatchesTransition.hh
 **
 ** \brief Top level component to control the analysis process.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_HASH_MATCH_FINDER_HH
#define iSAAC_ALIGNMENT_HASH_MATCH_FINDER_HH

#include <vector>

#include "alignment/Cluster.hh"
#include "alignment/Match.hh"
#include "reference/ReferenceHash.hh"

namespace isaac
{
namespace alignment
{

typedef reference::ContigList::Offset ReferenceOffset;
typedef std::vector<ReferenceOffset> ReferenceOffsetList;
typedef std::vector<ReferenceOffsetList> ReferenceOffsetLists;

template <typename ReferenceHash>
class SeedHashMatchFinder
{
protected:
    typedef typename ReferenceHash::KmerT KmerT;
    const ReferenceHash& referenceHash_;
    const unsigned seedBaseQualityMin_;

public:
    static const unsigned NO_CONTIG_FILTER = -1U;

    SeedHashMatchFinder(
        const ReferenceHash& referenceHash,
        const unsigned seedBaseQualityMin);

    template <bool reverse>
    std::size_t findSeedMatches(
        const Cluster& cluster,
        const flowcell::ReadMetadata& readMetadata,
        const std::size_t seedRepeatThreshold,
        Matches& matches,
        std::size_t& uncheckedSeeds) const;

    template <bool reverse, bool filterContigs>
    std::size_t getSeedAlignmentPositions(
        const reference::ContigList &contigList,
        const unsigned seedOffset,
        const flowcell::ReadMetadata &readMetadata,
        const unsigned filterContigId,
        typename ReferenceHash::MatchRange matchPositions,
        ReferenceOffsetList &referencePositions) const;
};

template <typename ReferenceHash, unsigned seedsPerMatchMax = 4>
class ClusterHashMatchFinder : SeedHashMatchFinder<ReferenceHash>
{
    typedef SeedHashMatchFinder<ReferenceHash> BaseT;
    typedef typename BaseT::KmerT KmerT;
public:
    using BaseT::NO_CONTIG_FILTER;
    static const unsigned SEED_LENGTH = ReferenceHash::SEED_LENGTH;
    static const unsigned SV_READ_SEEDS_MIN = 1;
    static const unsigned SHORT_READ_SEEDS_MIN = 2;
    static const unsigned LONG_READ_SEEDS_MIN = 3;
    BOOST_STATIC_ASSERT(seedsPerMatchMax >= SHORT_READ_SEEDS_MIN);

    ClusterHashMatchFinder(
        const ReferenceHash& referenceHash,
        const std::size_t candidateMatchesMax,
        const unsigned seedBaseQualityMin,
        const unsigned seedRepeatsMax);

    struct SeedHits
    {
        unsigned seedOffset_;
        typename ReferenceHash::MatchRange forwardMatches_;
        typename ReferenceHash::MatchRange reverseMatches_;

        std::size_t forwardHitCount() const
            {return std::distance(forwardMatches_.first, forwardMatches_.second);}
        std::size_t reverseHitCount() const
            {return std::distance(reverseMatches_.first, reverseMatches_.second);}
        std::size_t hitCount() const {return reverseHitCount() + forwardHitCount();}
        bool empty() const {return 0 == hitCount();}

        bool operator < (const SeedHits &that) const
        {
            return hitCount() < that.hitCount();
        }

        friend std::ostream &operator <<(std::ostream &os, const SeedHits &seedHits)
        {
            return os <<
                "SeedHits(" << seedHits.seedOffset_ << "," <<
                seedHits.forwardHitCount() << "," <<
                seedHits.reverseHitCount() << ")";
        }
    };

    std::size_t findReadMatches(
        const reference::ContigList &contigList,
        const Cluster& cluster,
        const flowcell::ReadMetadata& readMetadata,
        const std::size_t seedRepeatThreshold,
        MatchLists& matchLists,
        ReferenceOffsetLists& fwMergeBuffers,
        ReferenceOffsetLists& rvMergeBuffers) const;

    std::size_t findHeadAnchoredReadMatches(
        const reference::ContigList &contigList,
        const Cluster& cluster,
        const flowcell::ReadMetadata& readMetadata,
        const unsigned filterContigId,
        const std::size_t seedRepeatThreshold,
        const unsigned headSeedOffsetMax,
        MatchLists& matchLists,
        ReferenceOffsetLists& fwMergeBuffers,
        ReferenceOffsetLists& rvMergeBuffers) const;

private:
    const std::size_t candidateMatchesMax_;

    typedef common::StaticVector<SeedHits, 16> SeedsHits;

    template<bool filterContigs, bool detectStructuralVariant>
    void buildMatchesIteratively(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata& readMetadata,
        const Cluster& cluster,
        const unsigned filterContigId,
        const SeedsHits &seedsHits,
        MatchLists& matchLists,
        ReferenceOffsetLists& fwMergeBuffers,
        ReferenceOffsetLists& rvMergeBuffers) const;

    std::size_t collectSeedHits(
        const Cluster& cluster,
        const unsigned readIndex,
        const std::size_t seedRepeatThreshold,
        const unsigned endSeedOffset,
        SeedsHits& seedsHits) const;
};

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_HASH_MATCH_FINDER_HH
