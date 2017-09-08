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
 ** \file SeedGenerator.hh
 **
 ** \brief Generates kmers for a given permutation that match the mask
 **
 ** \author Roman Petrovski
 **/

#ifndef ISAAC_REFERENCE_SEED_GENERATOR_HH
#define ISAAC_REFERENCE_SEED_GENERATOR_HH

#include "oligo/KmerGenerator.hpp"
#include "reference/Contig.hh"
#include "reference/Seed.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
class SeedGeneratorThread: boost::noncopyable
{
public:
    /**
     * \param contigList index-ordered list of contigs
     */
    SeedGeneratorThread(const reference::ContigList &contigList) :
        contigList_(contigList)
    {
    }

protected:
    template <typename CallbackT>
    void thread(
        const unsigned threadNumber,
        const std::size_t threads,
        const CallbackT &callback) const;

private:
    /// index-ordered list of contigs
    const reference::ContigList &contigList_;
};


template <typename KmerT>
template <typename CallbackT>
void SeedGeneratorThread<KmerT>::thread(
    const unsigned threadNumber,
    const std::size_t threads,
    const CallbackT &callback) const
{
    typedef Seed<KmerT> SeedT;
    std::size_t generated = 0;
    for (const ContigList::Contig &contig : contigList_)
    {
        const std::size_t threadSectionLength = (contig.size() + threads - 1) / threads;

        const std::size_t beginGenomicOffset = threadSectionLength * threadNumber;
        if (contig.size() >= beginGenomicOffset + SeedT::SEED_LENGTH)
        {
            const std::size_t threadSectionLengthWithOverlap = threadSectionLength + SeedT::SEED_LENGTH - SeedT::STEP;
            const std::size_t threadSectionEnd = beginGenomicOffset + threadSectionLengthWithOverlap;
            const std::size_t endGenomicOffset = std::min(contig.size(), threadSectionEnd);

//            ISAAC_THREAD_CERR << "SeedGeneratorThread<KmerT>::thread " << threadNumber << " " << contig << " beginGenomicOffset: " << beginGenomicOffset << " endGenomicOffset: " << endGenomicOffset << std::endl;
            oligo::InterleavedKmerGenerator<SeedT::KMER_BASES, typename SeedT::KmerType, ContigList::Contig::const_iterator, SeedT::STEP> kmerGenerator(
                contig.begin() + beginGenomicOffset,
                contig.begin() + endGenomicOffset);

            typename SeedT::KmerType kmer(0);
            ContigList::Contig::const_iterator it;
            while (kmerGenerator.next(kmer, it))
            {
                const uint64_t kmerPosition = std::distance(contig.begin(), it);
//                ISAAC_THREAD_CERR << "kmerPosition:" << kmerPosition << std::endl;
                callback(threadNumber, kmer, contig.getIndex(), kmerPosition, false);
                ++generated;
            }
        }
    }
}

} // namespace reference
} // namespace isaac

#endif // #ifndef ISAAC_REFERENCE_PERMUTATED_KMER_GENERATOR_HH
