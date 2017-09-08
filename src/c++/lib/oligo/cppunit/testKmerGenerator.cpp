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
#include <boost/utility/binary.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testKmerGenerator.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestKmerGenerator, registryName("KmerGenerator"));

void TestKmerGenerator::setUp()
{
}

void TestKmerGenerator::tearDown()
{
}


void TestKmerGenerator::testSimple()
{
    typedef isaac::oligo::KmerGenerator<7, unsigned, std::vector<char>::const_iterator> KmerGenerator;
    unsigned kmer = 0;
    std::vector<char>::const_iterator position;
    {
        const std::string s("ANAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(kmer, 0x01B0U);
        CPPUNIT_ASSERT_EQUAL(position - v.begin(), 2L);
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAANAAAAAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AACGTAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(2L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAAAC"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(10L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAACG"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(11L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAACGT"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(12L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAACGTA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(13L, position - v.begin());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }

    {
        const std::string s = std::string("NAAACGTAAAAAAAAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAACGTA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(1L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AACGTAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(2L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("ACGTAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(3L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("CGTAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(4L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("GTAAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(5L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("TAAAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(6L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(7L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(8L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(9L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAAAC"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(10L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAAACG"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(11L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAACGT"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(12L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAACGTA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(13L, position - v.begin());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
}

void TestKmerGenerator::testStep2()
{
    typedef isaac::oligo::KmerGenerator<7, unsigned, std::vector<char>::const_iterator, 2> KmerGenerator;
    unsigned kmer = 0;
    std::vector<char>::const_iterator position;
    {
        const std::string s("ANAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAANAAAAAACGTGATTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin() + 1, v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAGGTA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(11L, position - v.begin());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }

}
void TestKmerGenerator::testInterleaved()
{
    typedef isaac::oligo::InterleavedKmerGenerator<7, unsigned, std::vector<char>::const_iterator, 2> KmerGenerator;
    unsigned kmer = 0;
    std::vector<char>::const_iterator position;
    {
        const std::string s("ANAACGTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }
    {
        const std::string s = std::string("ANAACGTAANAAAAAACGTGATTA");
        const std::vector<char> v(s.begin(), s.end());
        KmerGenerator kmerGenerator(v.begin(), v.end());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AACTAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(0L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("ACTAAAA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(2L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("CTAAAAC"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(4L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("TAAAACT"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(6L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAACTA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(8L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAACTAT"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(10L, position - v.begin());
        CPPUNIT_ASSERT(kmerGenerator.next(kmer, position));
        CPPUNIT_ASSERT_EQUAL(std::string("AAAGGTA"), std::string(isaac::oligo::bases<2, unsigned>(kmer, 7)));
        CPPUNIT_ASSERT_EQUAL(11L, position - v.begin());
        CPPUNIT_ASSERT(!kmerGenerator.next(kmer, position));
    }

}

