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
 ** \file ReferenceSorter.hh
 **
 ** Top level component to produce a sorted reference.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_HASH_HH
#define iSAAC_REFERENCE_REFERENCE_HASH_HH

#include <boost/format.hpp>
#include "common/NumaContainer.hh"
#include "oligo/Kmer.hh"

namespace isaac
{
namespace reference
{
template <typename KmerT> class ReferenceHasher;

template <typename KmerType, typename AllocatorT = std::allocator<void> >
class ReferenceHash
{
    typedef ReferenceHash<KmerType, AllocatorT> MyT;

public:
    typedef typename AllocatorT::template rebind<reference::ContigList::Offset> ReferenceOffsetAllocatorRebind;
    typedef typename ReferenceOffsetAllocatorRebind::other ReferenceOffsetAllocator;
    // for large genomes Offset type must be > 32 bit
    typedef reference::ContigList::Offset Offset;
    // offsets in linear genome indicating points where kmer is present
    typedef std::vector<Offset, ReferenceOffsetAllocator> Positions;
    typedef KmerType KmerT;
    static const unsigned SEED_LENGTH = oligo::KmerTraits<KmerT>::KMER_BASES;
    typedef typename Positions::const_iterator const_iterator;
    typedef std::pair<const_iterator, const_iterator> MatchRange;
    typedef void value_type;// compatibility with std containers for numa replications

    // the kmers are hashed into keys which are then used as indices into Offsets table
    typedef uint32_t KeyT;
    typedef typename AllocatorT::template rebind<Offset> OffsetAllocatorRebind;
    typedef typename OffsetAllocatorRebind::other OffsetAllocator;
    // offsets in Positions indicating ranges of offsets for the kmer
    // In reality if less than 4B kmers are unique, this one can get away with 32-bit
    // numbers even for genomes larger than 4B bases.
    typedef std::vector<Offset, OffsetAllocator> Offsets;

    KeyT keyFromKmer(KmerT kmer) const
    {
//        if (kmer == KmerT(/*0x02cee3cc14 */0x02e2c0c82b))
//        {
//            ISAAC_THREAD_CERR << "key:" << (kmer.bits_ % KEY_MAX) <<
//                "for kmer:" << oligo::Bases<oligo::BITS_PER_BASE, KmerT>(kmer, oligo::KmerTraits<KmerT>::KMER_BASES) << std::endl;
//        }

//        return kmer.bits_;
        return ((kmer.bits_ * a_ + b_) % largePrime_) % bucketCount_;
    }

    ReferenceHash(const uint64_t bucketCount)
        : a_(3308323), b_(7048005), largePrime_(1699023365707), bucketCount_(bucketCount), offsets_(bucketCount_, 0)
    {
        if (!bucketCount_)
        {
            BOOST_THROW_EXCEPTION(common::InvalidParameterException("Bucket count 0 is invalid"));
        }
        if (std::numeric_limits<KeyT>::max() < bucketCount_ - 1)
        {
            BOOST_THROW_EXCEPTION(common::InvalidParameterException(
                (boost::format("Bucket count %d is too large for key type %s. Max possible key value is %d") % bucketCount_ %
                    typeid(KeyT).name() % std::numeric_limits<KeyT>::max()
            ).str()));
        }
//        ISAAC_THREAD_CERR << "ReferenceHash()" << std::endl;
    }

    ReferenceHash(ReferenceHash &&that, const AllocatorT &allocator = AllocatorT())
        : a_(that.a_), b_(that.b_), largePrime_(that.largePrime_), bucketCount_(that.bucketCount_)
    {
        offsets_.swap(that.offsets_);
        positions_.swap(that.positions_);
//        ISAAC_THREAD_CERR << "ReferenceHash(ReferenceHash &&that, allocator)" << std::endl;
    }

    ReferenceHash(const ReferenceHash &that, const AllocatorT &allocator)
        : a_(that.a_), b_(that.b_), largePrime_(that.largePrime_), bucketCount_(that.bucketCount_)
        , offsets_(that.offsets_, allocator)
        , positions_(that.positions_, allocator)
    {
//        ISAAC_THREAD_CERR << "ReferenceHash(ReferenceHash &that, allocator)" << std::endl;
    }
//
//    void dumpDelta(const ReferenceHash &that)
//    {
//        for (KeyT k = 1; k < offsets_.size(); ++k)
//        {
//            if ((offsets_[k] - offsets_[k-1]) != (that.offsets_[k] - that.offsets_[k-1]))
//            {
//                int delta = int(offsets_[k] - offsets_[k-1]) - int(that.offsets_[k] - that.offsets_[k-1]);
//                if (abs(delta) > 1)
//                {
//                    ISAAC_THREAD_CERR << "delta>" << std::hex << std::setw(10) << k << " " << std::dec << int(offsets_[k] - offsets_[k-1]) - int(that.offsets_[k] - that.offsets_[k-1]) << "\n";
//                }
//            }
//        }
//    }
//
//    ReferenceHash &operator = (ReferenceHash &&that)
//    {
//        ISAAC_THREAD_CERR << "ReferenceHash &operator = (ReferenceHash &&that)" << std::endl;
//        a_ = that.a_;
//        b_ = that.b_;
//        largePrime_ = that.largePrime_;
//        bucketCount_ = that.bucketCount_;
//        offsets_.swap(that.offsets_);
//        positions_.swap(that.positions_);
//        return *this;
//    }

    MatchRange iSAAC_PROFILING_NOINLINE findMatches(const KmerT &kmer) const
    {
        const KeyT key = keyFromKmer(kmer);
        Offset positionsBegin = !key ? 0 : offsets_[key - 1];
        Offset positionsEnd = offsets_[key];
        ISAAC_ASSERT_MSG(positionsBegin <= positions_.size(), "Positions buffer overrun by positionsBegin:" << positionsBegin << " for kmer " << kmer);
        ISAAC_ASSERT_MSG(positionsBegin <= positionsEnd, "positionsEnd:" << positionsEnd << " overrun by positionsBegin:" << positionsBegin << " for kmer " << kmer);

        const MatchRange ret = std::make_pair(positions_.begin() + positionsBegin, positions_.begin() + positionsEnd);

    //    ISAAC_THREAD_CERR << "found " << std::distance(ret.first, ret.second) << " matches for " << oligo::Bases<oligo::BITS_PER_BASE, KmerT>(kmer, oligo::KmerTraits<KmerT>::KMER_BASES) << std::endl;
    //    BOOST_FOREACH(const ReferencePosition &pos, ret)
    //    {
    //        ISAAC_THREAD_CERR << pos << std::endl;
    //    }
        return ret;
    }

    MatchRange getEmptyRange() const
    {
        return std::make_pair(positions_.end(), positions_.end());
    }

    uint64_t getBucketCount() const {return bucketCount_;}
    uint64_t getA() const {return a_;}
    uint64_t getB() const {return b_;}
    uint64_t getLargePrime() const {return largePrime_;}
private:
    uint64_t a_;
    uint64_t b_;
    uint64_t largePrime_;
    uint64_t bucketCount_;
    Offsets offsets_;
//    std::vector<KmerT> uniqueKmers_;
    Positions positions_;

    friend class ReferenceHasher<MyT>;
};


template <typename HashType>
class NumaReferenceHash
{
    common::NumaContainerReplicas<HashType> replicas_;
public:
    typedef typename HashType::KmerT KmerT;
    typedef typename HashType::MatchRange MatchRange;
    typedef typename HashType::const_iterator const_iterator;
    typedef typename HashType::Positions Positions;
    typedef typename HashType::KeyT KeyT;
    typedef typename HashType::Offsets Offsets;
    static const unsigned SEED_LENGTH = HashType::SEED_LENGTH;

    NumaReferenceHash(HashType &&hash) :replicas_(std::move(hash))
    {
        ISAAC_THREAD_CERR << "NumaReferenceHash copy constructor" << std::endl;
    }

    MatchRange findMatches(const KmerT &kmer) const
    {
        return replicas_.threadNodeContainer().findMatches(kmer);
    }
};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_HASH_HH
