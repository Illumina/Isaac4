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

#ifndef iSAAC_REFERENCE_REFERENCE_HASHER_HH
#define iSAAC_REFERENCE_REFERENCE_HASHER_HH

#include "common/NumaContainer.hh"
#include "oligo/Kmer.hh"
//#include "PermutatedKmerGenerator.hh"
#include "SeedGenerator.hh"
#include "reference/ReferenceHash.hh"

namespace isaac
{
namespace reference
{

template <typename ReferenceHashT>
class ReferenceHasher
//    : PermutatedKmerGenerator<typename ReferenceHashT::KmerT, permutatedKmerGenerator::ForwardNoPermutate>
    : SeedGeneratorThread<typename ReferenceHashT::KmerT>
{
    typedef typename ReferenceHashT::KmerT KmerT;
    typedef SeedGeneratorThread<KmerT> BaseT;
//    typedef PermutatedKmerGenerator<KmerT, permutatedKmerGenerator::ForwardNoPermutate> BaseT;
    typedef typename ReferenceHashT::Positions Positions;
    typedef typename ReferenceHashT::Offset Offset;
    typedef typename ReferenceHashT::Offsets Offsets;
    static const std::size_t THREAD_BUFFER_KMERS_MAX = 8192; // arbitrary number that reduces the cost/benefit of acquiring a mutex
public:

    ReferenceHasher(const ContigList &contigList, common::ThreadVector &threads, const unsigned threadsMax);

    ReferenceHashT generate(const uint64_t bucketCount);
    void generate(ReferenceHashT &ret);

private:
    const ContigList &contigList_;
    common::ThreadVector &threads_;
    const unsigned threadsMax_;

    boost::ptr_vector<boost::mutex> mutexes_;
    typedef std::pair<typename ReferenceHashT::KmerT, ContigList::Offset> KmerWithPosition;
    typedef std::vector<KmerWithPosition> MutexBuffer;
    typedef std::vector<MutexBuffer> ThreadBuffer;
    typedef std::vector<ThreadBuffer> ThreadBuffers;
    ThreadBuffers threadBuffers_;

    void updateOffsets(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const KmerT &kmer);

    void countKmers(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const std::size_t threads);

    void storePosition(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const KmerT &kmer,
        const int contigId,
        const uint64_t kmerPosition,
        const bool reverse,
        const ContigList &contigList);

    void storePositions(
        ReferenceHashT &referenceHash,
        const unsigned threadNumber,
        const std::size_t threads,
        const ContigList &contigList);

    static void sortPositions(ReferenceHashT &referenceHash, const unsigned threadNumber, const std::size_t threads);
    static void updateEmptyOffsets(Offsets& offsets);
    static Offset countsToOffsets(Offsets& offsets);
    static void dumpCounts(boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash);
    static void dumpPositions(boost::mutex& mutex, MutexBuffer& buffer, ReferenceHashT& referenceHash);

    typename ReferenceHashT::KeyT mutexIdFromKey(const typename ReferenceHashT::KeyT key) const;
    void dumpDistribution(ReferenceHashT& ret);
};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_HASHER_HH
