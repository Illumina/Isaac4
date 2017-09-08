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
 ** \file Seed.hh
 **
 ** \brief Seeds are used to find alignment candidates in the reference genome
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_SEED_HH
#define iSAAC_REFERENCE_SEED_HH

#include <iostream>

#include "oligo/Kmer.hh"

namespace isaac
{
namespace reference
{

template <typename KmerT>
class InterleavedSeed
{
public:
    typedef KmerT KmerType;
    static const unsigned STEP = 1;
    static const unsigned KMER_BASES = KmerType::KMER_BASES;
    static const unsigned SEED_LENGTH = KMER_BASES * STEP;
    InterleavedSeed() : kmer_(0) {}
    InterleavedSeed(KmerT kmer) : kmer_(kmer) {}
    KmerType &kmer() {return kmer_;}
private:
    KmerType kmer_;
};


template <typename KmerT>
inline std::ostream &operator<<(std::ostream &os, const InterleavedSeed<KmerT> &seed)
{
    typedef InterleavedSeed<KmerT> SeedT;
    return os << "InterleavedSeed(" << oligo::Bases<oligo::BITS_PER_BASE, typename SeedT::KmerType>(seed.getKmer(), oligo::KmerTraits<typename SeedT::KmerType>::KMER_BASES) <<
        "(" << oligo::ReverseBases<oligo::BITS_PER_BASE, typename SeedT::KmerType>(seed.getKmer(), oligo::KmerTraits<typename SeedT::KmerType>::KMER_BASES) << ")" <<
         ")";
}

template <typename KmerT> struct Seed : public InterleavedSeed<KmerT>
{

};

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_SEED_HH
