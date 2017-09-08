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
 ** \file ReferenceHasher.cpp
 **
 ** Produces the mapping between kmer values and genomic positions
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "common/Exceptions.hh"
#include "common/Numa.hh"
#include "common/SystemCompatibility.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ReferenceHasher.hh"
#include "reference/SortedReferenceXml.hh"

namespace isaac
{
namespace reference
{

template <typename ReferenceHashT>
ReferenceHasher<ReferenceHashT>::ReferenceHasher (
    const ContigList &contigList,
    common::ThreadVector &threads,
    const unsigned threadsMax)
    : BaseT(contigList)
    , contigList_(contigList)
    , threads_(threads)
    , threadsMax_(threadsMax)
    , mutexes_(threadsMax_ * 2) // reduce collision probability somewhat
    , threadBuffers_(threadsMax_, ThreadBuffer(mutexes_.capacity()))

{
    while (mutexes_.capacity() != mutexes_.size())
    {
        mutexes_.push_back(new boost::mutex);
    }

    for (ThreadBuffer &threadBuffer : threadBuffers_)
    {
        for (MutexBuffer &mutexBuffer : threadBuffer)
        {
            mutexBuffer.reserve(THREAD_BUFFER_KMERS_MAX);
        }
    }
}

template<typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::dumpCounts(
    boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash)
{
    boost::lock_guard<boost::mutex> lock(mutex);
    for (const KmerWithPosition& kwp : buffer)
    {
        ++referenceHash.offsets_.at(referenceHash.keyFromKmer(kwp.first));
//        ++referenceHash.offsets_.at(kwp.first.bits_);
        //    ISAAC_THREAD_CERR << pos << std::endl;
    }
    buffer.clear();
}

template<typename ReferenceHashT>
typename ReferenceHashT::KeyT ReferenceHasher<ReferenceHashT>::mutexIdFromKey(const typename ReferenceHashT::KeyT key) const
{
    // we want the neighbor keys to be guarded by the same mutex.
    return (key / THREAD_BUFFER_KMERS_MAX) % mutexes_.size();
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::updateOffsets(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const KmerT &kmer)
{
    const typename ReferenceHashT::KeyT key = referenceHash.keyFromKmer(kmer);
    const std::size_t mutexIndex = mutexIdFromKey(key);
    MutexBuffer &buffer = threadBuffers_[threadNumber][mutexIndex];
    buffer.push_back(std::make_pair(kmer, ContigList::Offset(0)));

    if (buffer.capacity() == buffer.size())
    {
        dumpCounts(mutexes_.at(mutexIndex), buffer, referenceHash);
    }
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::countKmers(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const std::size_t threads)
{
    BaseT::thread(
        threadNumber, threads,
        [this, &referenceHash](
            const unsigned threadNumber, const KmerT &kmer, const unsigned contigIndex, const uint64_t kmerPosition, bool reverse)
        {
            updateOffsets(referenceHash, threadNumber, kmer);
        });

    ThreadBuffer &threadBuffer = threadBuffers_[threadNumber];
    for (std::size_t mutexIndex = 0; threadBuffer.size() > mutexIndex; ++mutexIndex)
    {
        dumpCounts(mutexes_.at(mutexIndex), threadBuffer[mutexIndex], referenceHash);
    }
}

template<typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::dumpPositions(
    boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash)
{
    boost::lock_guard<boost::mutex> lock(mutex);
    for (const KmerWithPosition& kwp : buffer)
    {
        const std::size_t offset = referenceHash.offsets_[referenceHash.keyFromKmer(kwp.first)]++;
        referenceHash.positions_.at(offset) = kwp.second;
//        ISAAC_THREAD_CERR << "kmer" << kwp.first << " offset:" << offset << " pos:" << kwp.second << std::endl;
    }
    buffer.clear();
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::storePosition(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const KmerT &kmer,
    const int contigId,
    const uint64_t kmerPosition,
    const bool reverse,
    const ContigList &contigList)
{
    const typename ReferenceHashT::KeyT key = referenceHash.keyFromKmer(kmer);
    const std::size_t mutexIndex = mutexIdFromKey(key);
    MutexBuffer &buffer = threadBuffers_[threadNumber][mutexIndex];
    ISAAC_ASSERT_MSG(!reverse, "This implementation does not support reverse kmers");
    buffer.push_back(std::make_pair(kmer, contigList.contigBeginOffset(contigId) + kmerPosition));

    if (buffer.capacity() == buffer.size())
    {
        dumpPositions(mutexes_.at(mutexIndex), buffer, referenceHash);
    }
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::storePositions(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const std::size_t threads,
    const ContigList &contigList)
{
    BaseT::thread(
        threadNumber, threads,
        [this, &referenceHash, &contigList](
            const unsigned threadNumber, const KmerT &kmer, const unsigned contigIndex, const uint64_t kmerPosition, bool reverse)
        {
            storePosition(referenceHash, threadNumber, kmer, contigIndex, kmerPosition, reverse, contigList);
        });

//    BaseT::kmerThread(
//        oligo::Permutate(oligo::KmerTraits<KmerT>::KMER_BASES), threadNumber, threads,
//        boost::bind(&ReferenceHasher::storePosition, this,
//                    boost::ref(referenceHash), _1, _2, _3, _4, _5,
//                    boost::ref(contigList)));
    ThreadBuffer &threadBuffer = threadBuffers_[threadNumber];
    for (std::size_t mutexIndex = 0; threadBuffer.size() > mutexIndex; ++mutexIndex)
    {
        dumpPositions(mutexes_.at(mutexIndex), threadBuffer[mutexIndex], referenceHash);
    }
}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::sortPositions(
    ReferenceHashT &referenceHash,
    const unsigned threadNumber,
    const std::size_t threads)
{
    const std::size_t blockLength = (referenceHash.offsets_.size() - 1) / threads + 1;
    const std::size_t blockBegin = blockLength * threadNumber;
    if (blockBegin < referenceHash.offsets_.size())
    {
        const std::size_t blockEnd = std::min(referenceHash.offsets_.size(), blockBegin + blockLength + 1);
        for (typename Offsets::iterator it = referenceHash.offsets_.begin() + blockBegin + 1;
            referenceHash.offsets_.begin() + blockEnd != it; ++it)
        {
            std::sort(referenceHash.positions_.begin() + *(it - 1), referenceHash.positions_.begin() + *it);
        }
    }
}

/**
 * \brief replaces 0 offsets with the offset of the previous k-mer
 */
template <typename KmerT>
void ReferenceHasher<KmerT>::updateEmptyOffsets(Offsets& offsets)
{
    Offset lastOffset = 0;
    BOOST_FOREACH(Offset &offset, offsets)
    {
        if (!offset)
        {
            offset = lastOffset;
        }
        else
        {
            lastOffset = offset;
        }
    }
}

/**
 * \brief replaces 0 offsets with the offset of the previous k-mer
 */
template <typename ReferenceHashT>
typename ReferenceHasher<ReferenceHashT>::Offset
ReferenceHasher<ReferenceHashT>::countsToOffsets(Offsets& offsets)
{
//    std::vector<Offset> top;
//    top.reserve(101);

    Offset offset = 0;
    for (Offset &count : offsets)
    {
        if (count)
        {
//            top.push_back(count);
//            std::push_heap(top.begin(), top.end(), [](const Offset &left, const Offset &right){return left > right;});
//            if (top.capacity() == top.size())
//            {
//                std::pop_heap(top.begin(), top.end(), [](const Offset &left, const Offset &right){return left > right;});
//                top.pop_back();
//            }
            using std::swap; swap(offset, count);
            offset += count;
        }
    }

//    std::sort(top.begin(), top.end(), [](const Offset &left, const Offset &right){return left < right;});
//    const Offset topTotal = std::accumulate(top.begin(), top.end(), 0);
//
////    for (const Offset &topOffset : top)
////    {
////        ISAAC_THREAD_CERR << "Top count:" << topOffset << " of " << topTotal << std::endl;
////    }
//    ISAAC_THREAD_CERR << "Top " << top.size() << " contain " << topTotal << " positions out of " << offset << std::endl;
    return offset;
}

//template <typename ReferenceHashT>
//std::size_t kmersToKeys(ReferenceHashT &referenceHash, typename ReferenceHashT::Offsets &offsets, std::vector<typename ReferenceHashT::KmerT> &uniqueKmers)
//{
//    const std::size_t uniqueCount = std::count_if(offsets.begin(), offsets.end(), [](const typename ReferenceHashT::KeyT &offset){return 0 != offset;});
//    uniqueKmers.clear();
//    uniqueKmers.reserve(uniqueCount);
//    std::size_t totalKmerCount = 0;
//    for (typename ReferenceHashT::KmerT kmer(0); offsets.size() > kmer.bits_; ++kmer.bits_)
//    {
//        if (offsets[kmer.bits_])
//        {
//            uniqueKmers.push_back(kmer);
//            totalKmerCount += offsets[kmer.bits_];
//        }
//    }
//
//    offsets.clear();
//    offsets.resize(referenceHash.getBucketCount(), 0);
//    for (typename ReferenceHashT::KmerT kmer : uniqueKmers)
//    {
//        ++offsets.at(referenceHash.keyFromKmer(kmer));
//    }
//    return totalKmerCount;
//}

template <typename ReferenceHashT>
ReferenceHashT ReferenceHasher<ReferenceHashT>::generate(const uint64_t bucketCount)
{
    ReferenceHashT ret(bucketCount);

    generate(ret);

    return ret;
}

//template<typename ReferenceHashT>
//void ReferenceHasher<ReferenceHashT>::dumpDistribution(ReferenceHashT& ret)
//{
//    //    const std::size_t uniqueKmers = std::count_if(ret.offsets_.begin(), ret.offsets_.end(), [](const typename ReferenceHashT::KeyT &offset){return 0 != offset;});
//    //    const std::size_t totalKmers = kmersToKeys(ret, ret.offsets_, ret.uniqueKmers_);
//    std::vector<unsigned> repeatDistribution(1000000, 0);
//    for (std::size_t key = 0; key < ret.offsets_.size(); ++key)
//    {
//        const unsigned count = ret.offsets_[key];
//        if (count)
//        {
//            //            ISAAC_THREAD_CERR << "Key>" << key << " " << count << "\n";
//            if (count >= repeatDistribution.size())
//            {
//                ++repeatDistribution[0];
//            }
//            else
//            {
//                ++repeatDistribution[count];
//            }
//        }
//    }
//    {
//        std::size_t total = 0;
//        for (unsigned repeatCount = 1; repeatDistribution.size() > repeatCount; ++repeatCount)
//        {
//            if (repeatDistribution[repeatCount])
//            {
//                total += repeatDistribution[repeatCount];
//                ISAAC_THREAD_CERR<< "RepeatDistribution>" << repeatCount << " " << repeatDistribution[repeatCount] << " " << total << std::endl;
//            }
//        }
//        total += repeatDistribution[0];
//        ISAAC_THREAD_CERR<< "RepeatDistribution>" << repeatDistribution.size() << " " << repeatDistribution[0] << " " << total << std::endl;
//    }
//}

template <typename ReferenceHashT>
void ReferenceHasher<ReferenceHashT>::generate(ReferenceHashT &ret)
{
    ISAAC_TRACE_STAT(
        "Constructing ReferenceHasher: for " << oligo::KmerTraits<KmerT>::KMER_BASES << "-mers ");

    threads_.execute(boost::bind(&ReferenceHasher::countKmers, this, boost::ref(ret), _1, _2), threadsMax_);

    static std::size_t maxUniqueKeys = 0;
//    const std::size_t uniqueKmers = std::count_if(ret.offsets_.begin(), ret.offsets_.end(), [](const typename ReferenceHashT::KeyT &offset){return 0 != offset;});
//    const std::size_t totalKmers = kmersToKeys(ret, ret.offsets_, ret.uniqueKmers_);

//    dumpDistribution(ret);
    const std::size_t uniqueKeys = std::count_if(ret.offsets_.begin(), ret.offsets_.end(), [](const typename ReferenceHashT::KeyT &offset){return 0 != offset;});
    maxUniqueKeys = std::max(maxUniqueKeys, uniqueKeys);
    const Offset total = countsToOffsets(ret.offsets_);
    ISAAC_THREAD_CERR <<
        " a:" << ret.getA() <<
        " b:" << ret.getB() <<
        " buckets:" << ret.getBucketCount() <<
        " and " << total <<
        " genome " << oligo::KmerTraits<KmerT>::KMER_BASES <<
        "-mers "
//        " and " << uniqueKmers <<
        " unique k-mers found " << uniqueKeys << " unique keys. maxUniqueKeys:" << maxUniqueKeys << std::endl;

//    return ret;
    ret.positions_.resize(total);
    ISAAC_TRACE_STAT(" reserving memory done for " << ret.positions_.size() << " positions");

    threads_.execute(
        [this, &ret](const unsigned threadNumber, const std::size_t threads)
        {
            ReferenceHasher::storePositions(ret, threadNumber, threads, contigList_);
        }, threadsMax_);

    ISAAC_THREAD_CERR << " generated " << total << " positions" << std::endl;

    updateEmptyOffsets(ret.offsets_);

    threads_.execute(
        [this, &ret](const unsigned threadNumber, const std::size_t threads)
        {
            sortPositions(ret, threadNumber, threads);
        }, threadsMax_);

    ISAAC_THREAD_CERR << " sorted " << ret.offsets_.back() << " positions" << std::endl;
}
//
template class ReferenceHasher<ReferenceHash<oligo::VeryShortKmerType> >;
//template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<16>, common::NumaAllocator<void, 0> > >;
//template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<17>, common::NumaAllocator<void, 0> > >;
//template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<18>, common::NumaAllocator<void, 0> > >;
//template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<19>, common::NumaAllocator<void, 0> > >;
//template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<20>, common::NumaAllocator<void, 0> > >;

template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<10>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<11>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<12>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<13>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<14>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<15>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<16>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<17>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<18>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<19>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;
template class ReferenceHasher<ReferenceHash<oligo::BasicKmerType<20>, common::NumaAllocator<void, common::numa::defaultNodeInterleave> > >;

} // namespace reference
} // namespace isaac
