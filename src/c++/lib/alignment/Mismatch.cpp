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
 ** \file Mismatch.cpp
 **
 ** No-intrinsics vectorization for mismatch counting.
 **
 ** \author Roman Petrovski
 **/

#include <stdint.h>
#include <string.h>

#include "common/SystemCompatibility.hh"

namespace isaac
{
namespace alignment
{

static const vint128_t o[] = {
    {0x0101010101010101, 0x0101010101010101},
    {0x0101010101010101 << 1, 0x0101010101010101 << 1},
    {0x0101010101010101 << 2, 0x0101010101010101 << 2},
    {0x0101010101010101 << 3, 0x0101010101010101 << 3},
    {0x0101010101010101 << 4, 0x0101010101010101 << 4},
    {0x0101010101010101 << 5, 0x0101010101010101 << 5},
    {0x0101010101010101 << 6, 0x0101010101010101 << 6},
    {uint64_t(0x0101010101010101) << 7, uint64_t(0x0101010101010101) << 7},
};

unsigned countMismatchesFast(
    const char* sequenceBegin,
    const char* sequenceEnd,
    const char* referenceBegin)
{
//    ISAAC_THREAD_CERR << "countMismatchesFast" << std::endl;
//
    unsigned ret = 0;
    while (sequenceBegin + sizeof(vint128_t) <= sequenceEnd)
    {
        size_t io = 0;
        vint128_t my = {0};
        while (sizeof(o) / sizeof(o[0]) > io && sequenceBegin + sizeof(vint128_t) <= sequenceEnd)
        {
            vint128_t  a;// = (const vint128_t *)&*sequenceBegin; // this leads to movdqa which will segfault
//            a = _mm_loadu_si128((const vint128_t *)&*sequenceBegin);
            memcpy((char*)&a, sequenceBegin, sizeof(a)); //_mm_loadu_si128((const vint128_t *)&*sequenceBegin);
            vint128_t  b;// = *(const vint128_t *)&*referenceBegin;
//            b = _mm_loadu_si128((const vint128_t *)&*referenceBegin);
            memcpy((char*)&b, referenceBegin, sizeof(b)); //_mm_loadu_si128((const vint128_t *)&*referenceBegin);

            vint128_t cmask;
//            cmask = _mm_cmpeq_epi8(a, b);
            for (unsigned i = 0; i < sizeof(vint128_t); ++i)
            {
                ((char*)&cmask)[i] = ((char*)&a)[i] == ((char*)&b)[i] ? 0xff:0x00;
            }

            vint128_t mask;
//            mask = _mm_andnot_si128(cmask, o[io]);
            mask = (~cmask & o[io]); //= _mm_andnot_si128(cmask, o[io]); // gcc 4.7 makes xor+and instead of andnot
//            my.m128i_ = _mm_or_si128(mask.m128i_, my.m128i_);
            my |= mask;//_mm_or_si128(mask.m128i_, my.m128i_);


            ++io;
            sequenceBegin += sizeof(vint128_t);
            referenceBegin += sizeof(vint128_t);
        }

        const uint64_t *masks = (uint64_t*)&my;
        ret +=__builtin_popcountll(*masks) + __builtin_popcountll(*(masks+1));

//        const unsigned int *uiMask = (unsigned int*)&my;
//        ret += __builtin_popcount(*uiMask) +
//            __builtin_popcount(*(uiMask + 1)) +
//            __builtin_popcount(*(uiMask + 2)) +
//            __builtin_popcount(*(uiMask + 3));

    }

    // After some testing turns out that vectorizing code above provides only 1-2% improvement on 2x100 data.
    // Just comment out the above if it cause trouble on a particular architecture

    // count the remainder of mismatches
    for (;sequenceEnd != sequenceBegin; ++sequenceBegin, ++referenceBegin)
    {
        ret += *sequenceBegin != *referenceBegin;
    }

    return ret;
}

} // namespace alignment
} // namespace isaac
