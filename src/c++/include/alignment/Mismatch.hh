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
 ** \file Mismatch.hh
 **
 ** \brief Basic alignment constants and utilities.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MISMATCH_HH
#define iSAAC_ALIGNMENT_MISMATCH_HH

#include <string>
#include <vector>
#include <stdint.h>

#include <boost/foreach.hpp>
#include <boost/ref.hpp>

#include "alignment/Read.hh"
#include "alignment/Quality.hh"
#include "common/FastIo.hh"
#include "oligo/Nucleotides.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{

/**
 * \brief defines match for the purpose of the alignment.
 */
inline bool isMatch(const char readBase, const char referenceBase)
{
//    return readBase == oligo::SEQUENCE_OLIGO_N || (readBase == referenceBase && referenceBase != oligo::REFERENCE_OLIGO_N);

    // reference N should not match (unless corresponding sequence is N as well) because then stretches of N in reference attract
    // piles of reads to positions where reference begins to have some real bases.
    // sequence N should match anything because variable length input reads are padded with (sometimes) many N. If we count those
    // N as mismatches, all padded reads will have bad MAPQs.
    // TODO: Check what proportion of N is needed for this to happen.
    // TODO: Check if the assumption is valid at all as N qscores are very low anyway.
//    return readBase == oligo::SEQUENCE_OLIGO_N || readBase == referenceBase;
    return readBase == referenceBase;
}

inline bool isMismatch(const char readBase, const char referenceBase)
{
    return !isMatch(readBase, referenceBase);
}

unsigned countMismatchesFast(
    const char* sequenceBegin,
    const char* sequenceEnd,
    const char* referenceBegin);

inline unsigned iSAAC_PROFILING_NOINLINE countMismatches(
    std::vector<char>::const_iterator sequenceBegin,
    std::vector<char>::const_iterator sequenceEnd,
    reference::Contig::const_iterator referenceBegin)
{
    const unsigned fastMismatches = countMismatchesFast(&*sequenceBegin, &*sequenceEnd, &*referenceBegin);
    return fastMismatches;
//    const unsigned mismatches =
//        std::inner_product(sequenceBegin, sequenceEnd, referenceBegin, 0, std::plus<unsigned>(), &isMismatch);
//    return mismatches;
}


/**
 * \brief moves sequenceBegin to the first position followed by CONSECUTIVE_MATCHES_MAX matches
 *
 * \return pair(distance moved, edit distance adjustment). edit distance adjustment equals all mismatches that have been clipped away.
 *         Note that N is considered to be an edit distance mismatch in this case.
 */
template <unsigned CONSECUTIVE_MATCHES_MIN, typename SequenceIteratorT, typename ReferenceIteratorT, typename BaseExtractor>
std::pair<unsigned, unsigned> clipMismatches(
    SequenceIteratorT sequenceBegin, const SequenceIteratorT sequenceEnd,
    ReferenceIteratorT referenceBegin, ReferenceIteratorT referenceEnd,
    BaseExtractor baseExtractor)
{
    unsigned matchesInARow = 0;
    unsigned ediDistanceMismatches = 0;
    // The number of mismatches that are not part of the sequence that gets clipped
    unsigned ediDistanceMismatchesUnclipped = 0;
    unsigned ret = 0;
    while (sequenceEnd != sequenceBegin && referenceBegin != referenceEnd && CONSECUTIVE_MATCHES_MIN > matchesInARow)
    {
        char sequenceBase = baseExtractor(*sequenceBegin);
        if (isMatch(sequenceBase, *referenceBegin))
        {
            ++matchesInARow;
            ediDistanceMismatchesUnclipped += (sequenceBase != *referenceBegin);
        }
        else
        {
            matchesInARow = 0;
            ediDistanceMismatchesUnclipped = 0;
        }
        ediDistanceMismatches += (sequenceBase != *referenceBegin);
        ++sequenceBegin;
        ++referenceBegin;
        ++ret;
    }

    return (CONSECUTIVE_MATCHES_MIN == matchesInARow) ?
        std::make_pair(ret - matchesInARow, ediDistanceMismatches - ediDistanceMismatchesUnclipped):
        std::make_pair(0U,0U);
}

template <typename SequenceIteratorT, typename BaseExtractor>
unsigned countMatches(
    SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    reference::Contig::const_iterator referenceBegin,
    const reference::Contig::const_iterator referenceEnd,
    BaseExtractor baseExtractor)
{
    unsigned ret = 0;
    for (;sequenceEnd != sequenceBegin && referenceEnd != referenceBegin;
        ++sequenceBegin, ++referenceBegin)
    {
        ret += isMatch(baseExtractor(*sequenceBegin), *referenceBegin);
    }
    return ret;
}

template <typename SequenceIteratorT>
unsigned countMatches(
    const SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    const reference::Contig::const_iterator referenceBegin,
    const reference::Contig::const_iterator referenceEnd)
{
    return countMatches(sequenceBegin, sequenceEnd,
                        referenceBegin, referenceEnd,
                        [](typename std::iterator_traits<SequenceIteratorT>::value_type c){return c;});
}

template <typename SequenceIteratorT, typename ReferenceIteratorT, typename BaseExtractor>
unsigned countMismatches(
    SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    ReferenceIteratorT referenceBegin,
    const ReferenceIteratorT referenceEnd,
    BaseExtractor baseExtractor)
{
//    ISAAC_ASSERT_MSG(std::distance(sequenceBegin, sequenceEnd) < 1000, "tada " << std::distance(sequenceBegin, sequenceEnd));
//    ISAAC_THREAD_CERR << " countMismatches " << std::distance(sequenceBegin, sequenceEnd) << "compareLength " <<
//        " read '" << common::makeFastIoString(sequenceBegin, sequenceEnd) <<
//        "' ref '" << common::makeFastIoString(referenceBegin, referenceBegin + std::distance(sequenceBegin, sequenceEnd)) << "'" << std::endl;

    unsigned ret = 0;
    for (;sequenceEnd != sequenceBegin && referenceEnd != referenceBegin;
        ++sequenceBegin, ++referenceBegin)
    {
        ret += !isMatch(baseExtractor(*sequenceBegin), *referenceBegin);
    }

    return ret;
}

template <typename SequenceIteratorT, typename ReferenceIteratorT, typename BaseExtractor>
unsigned firstMismatchOffset(
    const SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    ReferenceIteratorT referenceBegin,
    const ReferenceIteratorT referenceEnd,
    BaseExtractor baseExtractor)
{
//    ISAAC_ASSERT_MSG(std::distance(sequenceBegin, sequenceEnd) < 1000, "tada " << std::distance(sequenceBegin, sequenceEnd));
//    std::cerr << " firstMismatchOffset " << std::distance(sequenceBegin, sequenceEnd) <<
//        " read '" << common::makeFastIoString(sequenceBegin, sequenceEnd) <<
//        "' ref '" << common::makeFastIoString(referenceBegin, referenceBegin + std::distance(sequenceBegin, sequenceEnd)) << "'" << std::endl;

    SequenceIteratorT sequence = sequenceBegin;
    for (;sequenceEnd != sequence && referenceEnd != referenceBegin;
        ++sequence, ++referenceBegin)
    {
        if (!isMatch(baseExtractor(*sequence), *referenceBegin))
        {
            break;
        }
    }

    return std::distance(sequenceBegin, sequence);
}



template <typename SequenceIteratorT, typename ReferenceIteratorT>
unsigned countMismatches(
    const SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd,
    const ReferenceIteratorT referenceBegin,
    const ReferenceIteratorT referenceEnd)
{
    return countMismatches(sequenceBegin, sequenceEnd,
                           referenceBegin, referenceEnd,
                           [](typename std::iterator_traits<SequenceIteratorT>::value_type c){return c;});
}


template <typename SequenceIteratorT, typename ReferenceIteratorT, typename BaseExtractor>
unsigned countMismatches(
    const SequenceIteratorT basesIterator,
    const ReferenceIteratorT referenceBegin,
    const ReferenceIteratorT referenceEnd,
    int length,
    BaseExtractor baseExtractor)
{
    ISAAC_ASSERT_MSG(0 <= length, "Positive length is required:" << length);
    return countMismatches(basesIterator, basesIterator + length,
                           referenceBegin, referenceEnd, baseExtractor);
}

/**
 * \brief counts the number of mismatches between the reference and sequence. Unlike the ones above,
 *        consider any discrepancy between the reference and sequence to be a mismatch
 */
template <typename SequenceIteratorT>
inline unsigned countEditDistanceMismatches(
    const reference::ContigList &reference,
    const SequenceIteratorT basesIterator,
    const reference::ReferencePosition pos,
    unsigned length)
{
    const reference::Contig &contig = reference.at(pos.getContigId());
    reference::Contig::const_iterator referenceBaseIt = contig.begin() + pos.getPosition();
    const unsigned compareLength = std::min<unsigned>(length, std::distance(referenceBaseIt, contig.end()));
    unsigned mismatches = 0;
    BOOST_FOREACH(const unsigned char readBase, std::make_pair(basesIterator, basesIterator + compareLength))
    {
        mismatches += *referenceBaseIt != oligo::getReferenceBaseFromBcl(readBase);
        ++referenceBaseIt;
    }

    referenceBaseIt = contig.begin() + pos.getPosition();

//    ISAAC_THREAD_CERR << mismatches << " mismatches " << compareLength << "compareLength " << pos <<
//        " read '" << oligo::bclToString(basesIterator, compareLength) <<
//        "' ref '" << common::makeFastIoString(referenceBaseIt, referenceBaseIt + compareLength) << "'" <<
//        std::endl;

    return mismatches;
}

template <typename SequenceIteratorT, std::size_t HOMOPOLYMER_LENGTH_MIN = 16>
bool containsHomopolymer(
    SequenceIteratorT sequenceBegin,
    const SequenceIteratorT sequenceEnd)
{
    char lastBase = 0;
    std::size_t repeats = 0;
    while (sequenceBegin != sequenceEnd)
    {
        if (*sequenceBegin == lastBase)
        {
            ++repeats;
            if (HOMOPOLYMER_LENGTH_MIN == repeats)
            {
                return true;
            }
        }
        else
        {
            lastBase = *sequenceBegin;
            repeats = 0;
        }
        ++sequenceBegin;
    }
    return false;
}

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MISMATCH_HH
