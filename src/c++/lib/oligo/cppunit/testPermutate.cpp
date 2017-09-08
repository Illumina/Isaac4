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
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/math/tools/big_constant.hpp>
using namespace std;

#include "RegistryName.hh"
#include "testPermutate.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestPermutate, registryName("Permutate"));

namespace isaac
{
namespace oligo
{
    template <> struct KmerBitsType<48>{typedef __uint128_t BitsType;};
}
}

void TestPermutate::setUp()
{
}

void TestPermutate::tearDown()
{
}

void TestPermutate::testFourBlocks()
{
    using namespace isaac;
    using boost::assign::list_of;
    const oligo::KmerType kmer(0xFEDCBA9876543210UL);
    const oligo::KmerType abcd(0xFEDCBA9876543210UL);
    const oligo::KmerType adbc(0xFEDC3210BA987654UL);
    const oligo::KmerType dbca(0x3210BA987654FEDCUL);
    const unsigned blockLength = 8;
    const std::vector<unsigned char> ABCD = list_of(0)(1)(2)(3);
    const std::vector<unsigned char> ADBC = list_of(0)(3)(1)(2);
    const std::vector<unsigned char> DBCA = list_of(3)(1)(2)(0);
    const oligo::Permutate ABCD_ABCD(blockLength, ABCD, ABCD, 0);
    const oligo::Permutate ABCD_ADBC(blockLength, ABCD, ADBC, 0);
    const oligo::Permutate ABCD_DBCA(blockLength, ABCD, DBCA, 0);
    const oligo::Permutate ADBC_DBCA(blockLength, ADBC, DBCA, 0);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(kmer), oligo::KmerType(0xFEDC3210BA987654UL));
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(kmer), oligo::KmerType(0x3210BA987654FEDCUL));
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(kmer), oligo::KmerType(0xBA9876543210FEDCUL));
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD(abcd), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC(abcd), adbc);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA(abcd), dbca);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA(adbc), dbca);
    CPPUNIT_ASSERT_EQUAL(ABCD_ABCD.reorder(abcd), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_ADBC.reorder(adbc), abcd);
    CPPUNIT_ASSERT_EQUAL(ABCD_DBCA.reorder(dbca), abcd);
    CPPUNIT_ASSERT_EQUAL(ADBC_DBCA.reorder(dbca), abcd);
}

void TestPermutate::testEightBlocks()
{
    using namespace isaac;
    using boost::assign::list_of;
    // A  B  C  D  E  F  G  H
    //FE DC BA 98 76 54 32 10
    const oligo::KmerType kmer = oligo::KmerType(0xFEDCBA9876543210UL);
    const oligo::KmerType abcdefgh = oligo::KmerType(0xFEDCBA9876543210UL);
    const oligo::KmerType adbcefgh = oligo::KmerType(0xFE98DCBA76543210UL);
    const oligo::KmerType dghbcafe = oligo::KmerType(0x983210DCBAFE5476UL);
    const unsigned blockLength = 4;
    const std::vector<unsigned char> ABCDEFGH = list_of(0)(1)(2)(3)(4)(5)(6)(7);
    const std::vector<unsigned char> ADBCEFGH = list_of(0)(3)(1)(2)(4)(5)(6)(7);
    const std::vector<unsigned char> DGHBCAFE = list_of(3)(6)(7)(1)(2)(0)(5)(4);
    const oligo::Permutate ABCDEFGH_ABCDEFGH(blockLength, ABCDEFGH, ABCDEFGH, 0);
    const oligo::Permutate ABCDEFGH_ADBCEFGH(blockLength, ABCDEFGH, ADBCEFGH, 0);
    const oligo::Permutate ABCDEFGH_DGHBCAFE(blockLength, ABCDEFGH, DGHBCAFE, 0);
    const oligo::Permutate ADBCEFGH_DGHBCAFE(blockLength, ADBCEFGH, DGHBCAFE, 0);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(kmer), kmer);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(kmer), oligo::KmerType(0xFE98DCBA76543210UL));
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(kmer), oligo::KmerType(0x983210DCBAFE5476UL));
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(kmer), oligo::KmerType(0xDC3210BA98FE5476UL));
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH(abcdefgh), adbcefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE(abcdefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE(adbcefgh), dghbcafe);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ABCDEFGH.reorder(abcdefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_ADBCEFGH.reorder(adbcefgh), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ABCDEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
    CPPUNIT_ASSERT_EQUAL(ADBCEFGH_DGHBCAFE.reorder(dghbcafe), abcdefgh);
}

template <typename KmerT>
void testPermutate(KmerT original, const KmerT expected, const std::vector<isaac::oligo::Permutate> permutateList)
{
    CPPUNIT_ASSERT_EQUAL(original, permutateList.front()(original));
    CPPUNIT_ASSERT_EQUAL(original, permutateList.front().reorder(original));
    KmerT permuted = original;
    BOOST_FOREACH(const isaac::oligo::Permutate &permutate, permutateList)
    {
        permuted = permutate(permuted);
        CPPUNIT_ASSERT_EQUAL(original, permutate.reorder(permuted));
    }
    CPPUNIT_ASSERT_EQUAL(expected, permuted);
    CPPUNIT_ASSERT_EQUAL(original, permutateList.back().reorder(permuted));
}


void TestPermutate::testNoErrors()
{
    const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
    const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);

    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(0, true, false);
        CPPUNIT_ASSERT_EQUAL(1UL, permutateList.size());
        testPermutate(ORIGINAL16, ORIGINAL16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(0, true, false);
        CPPUNIT_ASSERT_EQUAL(1UL, permutateList.size());
        testPermutate(ORIGINAL, ORIGINAL, permutateList);
    }
}


void TestPermutate::testTwoErrors()
{
    const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
    const isaac::oligo::ShortKmerType EXPECTED16(0x10547632U);
    const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);
    const isaac::oligo::KmerType EXPECTED(0x3210ba98fedc7654UL);

    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(2, 0, false);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL16, EXPECTED16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(2, 0, false);
        CPPUNIT_ASSERT_EQUAL(6UL, permutateList.size());
        testPermutate(ORIGINAL, EXPECTED, permutateList);
    }

}

void TestPermutate::testThreeErrors()
{
    const isaac::oligo::ShortKmerType ABCDEFGH16(0x76543210U);
    const isaac::oligo::ShortKmerType FGHDEABC16(0x01543762U);
    const isaac::oligo::KmerType ABCDEFGH(0xFEDCBA9876543210UL);
    const isaac::oligo::KmerType FGHDEABC(0x1032ba9876fedc54UL);

    using namespace isaac;
    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(3, 0, false);
        CPPUNIT_ASSERT_EQUAL(56UL, permutateList.size());
        testPermutate(ABCDEFGH16, FGHDEABC16, permutateList);
    }

    {
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(3, 0, false);
        CPPUNIT_ASSERT_EQUAL(56UL, permutateList.size());
        testPermutate(ABCDEFGH, FGHDEABC, permutateList);
    }

}


void TestPermutate::testFourErrors()
{
    using namespace isaac;
    {
        const isaac::oligo::ShortKmerType ORIGINAL16(0x76543210U);
        const isaac::oligo::ShortKmerType EXPECTED16(0x02147653U);
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::ShortKmerType>(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL16, EXPECTED16, permutateList);
    }

    {
        const isaac::oligo::KmerType ORIGINAL(0xFEDCBA9876543210UL);
        const isaac::oligo::KmerType EXPECTED(0x10543298fedcba76UL);
        const std::vector<oligo::Permutate> permutateList = oligo::getPermutateList<oligo::KmerType>(4, 0, false);
        CPPUNIT_ASSERT_EQUAL(70UL, permutateList.size());
        testPermutate(ORIGINAL, EXPECTED, permutateList);
    }
}
