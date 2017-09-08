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
#include <algorithm>
#include <cstdlib>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp> 
#include <boost/lambda/construct.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testFragmentBuilder.hh"
#include "BuilderInit.hh"

#include "flowcell/SequencingAdapterMetadata.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestFragmentBuilder, registryName("FragmentBuilder"));


TestFragmentBuilder::TestFragmentBuilder()
    : readMetadataList(getReadMetadataList())
    , flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, isaac::flowcell::FastqFlowcellData(false, '!', false), 8, 0, std::vector<unsigned>(),
                                           readMetadataList, "blah"))

    , contigList(getContigList())
    , bcl0(getBcl(readMetadataList, contigList, 0, 2, 3))
    , bcl2(getBcl(readMetadataList, contigList, 2, 1, 2))
    , bcl3(getBclClusters(readMetadataList,
           subv(getBclVector(readMetadataList, contigList, 0, 2, 3), 0,4) +
           "x" +  // x == 30*4+0 -- replaces 'T' with 'A'
           subv(getBclVector(readMetadataList, contigList, 0, 2, 3), 5, 190) +
           "Q" +  // x == 20*4+0 -- replaces 'A' with 'C'
           subv(getBclVector(readMetadataList, contigList, 0, 2, 3), 196)))
    , bcl4l(getBclClusters(readMetadataList, getBcl(substr(contigList[4], 0, 44) + substr(contigList[4], 0, 56) +
                   substr(reverseComplement(contigList[4]), 0, 42) + substr(reverseComplement(contigList[4]), 0, 58))))
    , bcl4t(getBclClusters(readMetadataList,getBcl(substr(contigList[4], 16, 44) + substr(contigList[4], 0, 56) +
                   substr(reverseComplement(contigList[4]), 18, 42) + substr(reverseComplement(contigList[4]), 0, 58))))
    , bcl4lt(getBclClusters(readMetadataList,getBcl(std::vector<char>(10, 'A') + contigList[4] + std::vector<char>(30, 'C') +
                    std::vector<char>(15, 'G') + reverseComplement(contigList[4]) + std::vector<char>(25, 'T'))))
    , tile0(32)
    , tile2(31)
    , clusterId0(1234)
    , clusterId2(12345)
    , cluster0(getMaxReadLength(readMetadataList))
    , cluster2(getMaxReadLength(readMetadataList))
    , cluster3(getMaxReadLength(readMetadataList))
    , cluster4l(getMaxReadLength(readMetadataList))
    , cluster4t(getMaxReadLength(readMetadataList))
    , cluster4lt(getMaxReadLength(readMetadataList))
{
    cluster0.init(readMetadataList, bcl0.cluster(0), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0, 0);
    cluster2.init(readMetadataList, bcl2.cluster(0), tile2, clusterId2, isaac::alignment::ClusterXy(0,0), true, 0, 0);
    cluster3.init(readMetadataList, bcl3.cluster(0), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0, 0);
    cluster4l.init(readMetadataList, bcl4l.cluster(0), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0, 0);
    cluster4t.init(readMetadataList, bcl4t.cluster(0), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0, 0);
    cluster4lt.init(readMetadataList, bcl4lt.cluster(0), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0, 0);
}

void TestFragmentBuilder::setUp()
{
    matchList.clear();
}

void TestFragmentBuilder::tearDown()
{
}

//static const isaac::alignment::SequencingAdapter testAdapters;

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

void TestFragmentBuilder::testEmptyMatchList()
{
    using isaac::alignment::templateBuilder::FragmentBuilder;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    isaac::alignment::FragmentMetadataList fragments;
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, false, alignmentCfg, cigarBuffer, false);
    CPPUNIT_ASSERT(fragments.empty());
    CPPUNIT_ASSERT(cigarBuffer.empty());
}

void TestFragmentBuilder::auxSingleSeed(const unsigned s0, const unsigned s1)
{
    const unsigned offset0 = 0;
    const unsigned offset1 = 32;
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(2 + offset0, false));
    matchList.push_back(Match(175 - offset1, true));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)2, matchList.size());
//    CPPUNIT_ASSERT_EQUAL((uint64_t)s0, matchList[0].seedId.unusedgetSeed());
//    CPPUNIT_ASSERT_EQUAL((uint64_t)s1, matchList[1].seedId.unusedgetSeed());
//    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[0].seedId.unusedgetSeed()].getReadIndex());
//    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[1].seedId.unusedgetSeed()].getReadIndex());
    // Create the fragment builder
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 456, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 1, cluster0, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 1, matchList.end(), cluster0, true, fragments[1]);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragments.size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, cigarBuffer.size());
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragments[0].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)2, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[0][0].logProbability, (double)0.000001);
    // Check that the fragments actually align at the expected positions
    for (size_t i = 0; 100 > i; ++i)
    {
        const unsigned position = fragments[0][0].position;
        CPPUNIT_ASSERT_EQUAL(cluster0[0].getForwardSequence()[i], contigList[0][i + position]);
    }
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragments[1].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)107, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[1][0].logProbability, (double)0.000001);
    // Check that the fragments actually align at the expected positions
    const std::vector<char> reverse0 = reverseComplement(contigList[0]);
    for (size_t i = 0; 100 > i; ++i)
    {
        CPPUNIT_ASSERT_EQUAL(cluster0[1].getForwardSequence()[i], reverse0[i + 3]);
        const unsigned position = fragments[1][0].position;
        const unsigned length = fragments[1][0].getObservedLength();
        const char reverseBase = isaac::oligo::reverseBase(cluster0[1].getForwardSequence()[i]);
        CPPUNIT_ASSERT_EQUAL(reverseBase, contigList[0][position + length - 1 - i]);
    }
}

void TestFragmentBuilder::testSingleSeed()
{
    // test on seed index 0 and 3 (both at offset 0)
    auxSingleSeed(0, 3);
}

void TestFragmentBuilder::testSeedOffset()
{
    // test on seed index 0 and 3 (both at offset 0)
    auxSingleSeed(1, 5);
}

void TestFragmentBuilder::testMultiSeed()
{
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(2, false));
    matchList.push_back(Match(2 + 32, false));
    matchList.push_back(Match(2 + 64, false));
    matchList.push_back(Match(175, true));
    matchList.push_back(Match(175 - 32, true));
    matchList.push_back(Match(175 - 64, true));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)6, matchList.size());
    // Create the fragment builder
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 3, cluster0, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 3, matchList.end(), cluster0, true, fragments[1]);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragments.size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, cigarBuffer.size());
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragments[0].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)2, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[0][0].logProbability, (double)0.000001);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragments[1].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)107, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[1][0].logProbability, (double)0.000001);
}

void TestFragmentBuilder::testRepeats()
{
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;

    matchList.push_back(Match(contigList.contigBeginOffset(2) + 1, false));
    matchList.push_back(Match(contigList.contigBeginOffset(3) + 6, false));
    matchList.push_back(Match(contigList.contigBeginOffset(2) + 1 + 32, false));
    matchList.push_back(Match(contigList.contigBeginOffset(2) + 1 + 64, false));
    matchList.push_back(Match(contigList.contigBeginOffset(3) + 6 + 64, false));
    matchList.push_back(Match(contigList.contigBeginOffset(3) + 201, true));
    matchList.push_back(Match(contigList.contigBeginOffset(2) + 196, true));
    matchList.push_back(Match(contigList.contigBeginOffset(2) + 196 - 32, true));
    matchList.push_back(Match(contigList.contigBeginOffset(3) + 201 - 32, true));
    matchList.push_back(Match(contigList.contigBeginOffset(2) + 196 - 64, true));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)10, matchList.size());
    // Create the fragment builder
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 5, cluster2, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 5, matchList.end(), cluster2, true, fragments[1]);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragments.size());
    CPPUNIT_ASSERT_EQUAL((size_t)4, cigarBuffer.size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragments[0].size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragments[1].size());
    // First fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)1, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[0][0].logProbability, (double)0.000001);
    // Second fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragments[0][1].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)6, fragments[0][1].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[0][1].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][1].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][1].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][1].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][1].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][1].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[0][1].logProbability, (double)0.000001);
    // First fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)128, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[1][0].logProbability, (double)0.000001);
    // Seconf fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragments[1][1].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)133, fragments[1][1].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[1][1].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][1].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][1].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragments[1][1].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][1].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][1].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-0.0100005, fragments[1][1].logProbability, (double)0.000001);
}

void TestFragmentBuilder::testMismatches()
{
    CPPUNIT_ASSERT_EQUAL(bcl0.getClusterCount() * bcl0.getClusterLength(), bcl3.getClusterCount() * bcl3.getClusterLength());
    // Check the mismatch on the forward strand
    CPPUNIT_ASSERT_EQUAL(*(bcl0.cluster(0) + 3), *(bcl3.cluster(0) + 3));
    CPPUNIT_ASSERT((*(bcl0.cluster(0) + 4)&3) != (*(bcl3.cluster(0) + 4)&3));
    CPPUNIT_ASSERT_EQUAL(*(bcl0.cluster(0) + 5), *(bcl3.cluster(0) + 5));
    // Check the mismatch on the reverse strand
    CPPUNIT_ASSERT_EQUAL(*(bcl0.cluster(0) + 1944), *(bcl3.cluster(0) + 194));
    CPPUNIT_ASSERT((*(bcl0.cluster(0) + 195)&3) != (*(bcl3.cluster(0) + 195)&3));
    CPPUNIT_ASSERT_EQUAL(*(bcl0.cluster(0) + 196), *(bcl3.cluster(0) + 196));
//    const unsigned s0 = 0;
//    const unsigned s1 = 3;
    const unsigned offset0 = 0;
    const unsigned offset1 = 32;
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(contigList.contigBeginOffset(0) + 2 + offset0, false));
    matchList.push_back(Match(contigList.contigBeginOffset(0) + 175 - offset1, true));
    // Sanity check on the match list
    CPPUNIT_ASSERT_EQUAL((size_t)2, matchList.size());
//    CPPUNIT_ASSERT_EQUAL((uint64_t)s0, matchList[0].seedId.unusedgetSeed());
//    CPPUNIT_ASSERT_EQUAL((uint64_t)s1, matchList[1].seedId.unusedgetSeed());
//    CPPUNIT_ASSERT_EQUAL((unsigned)0, seedMetadataList[matchList[0].seedId.unusedgetSeed()].getReadIndex());
//    CPPUNIT_ASSERT_EQUAL((unsigned)1, seedMetadataList[matchList[1].seedId.unusedgetSeed()].getReadIndex());
    // Create the fragment builder
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 1, cluster3, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 1, matchList.end(), cluster3, true, fragments[1]);
    // check buffer geometry
    CPPUNIT_ASSERT_EQUAL((size_t)2, fragments.size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, cigarBuffer.size());
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragments[0].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)2, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[0][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-8.016268063, fragments[0][0].logProbability, (double)0.000000001);
    // Check that the fragments actually align at the expected positions
    CPPUNIT_ASSERT_EQUAL(cluster3[0].getForwardSequence()[4], 'A');
    for (size_t i = 0; 100 > i; ++i)
    {
        if (4 != i)
        {
            CPPUNIT_ASSERT_EQUAL(cluster3[0].getForwardSequence()[i], contigList[0][i + 2]);
        }
    }
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((size_t)1, fragments[1].size());
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].contigId);
    CPPUNIT_ASSERT_EQUAL((int64_t)107, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(100U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)1, fragments[1][0].mismatchCount);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)-5.713682970, fragments[1][0].logProbability, (double)0.000000001);
    // Check that the fragments actually align at the expected positions
    const std::vector<char> reverse0 = reverseComplement(contigList[0]);
    CPPUNIT_ASSERT(cluster3[1].getForwardSequence()[95] != reverse0[95 + 3]);
    for (size_t i = 0; 100 > i; ++i)
    {
        if (95 != i)
        {
            CPPUNIT_ASSERT_EQUAL(cluster3[1].getForwardSequence()[i], reverse0[i + 3]);
        }
    }
}

void TestFragmentBuilder::testLeadingSoftClips()
{
//    const unsigned s0 = 2;
//    const unsigned s1 = 5;
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(contigList.contigBeginOffset(4) + 20, false));
    matchList.push_back(Match(contigList.contigBeginOffset(4) + 6, true));
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 1, cluster4l, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 1, matchList.end(), cluster4l, true, fragments[1]);
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((int64_t)0, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(56U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((44<<4)|Cigar::SOFT_CLIP), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((56<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)44, fragments[0][0].mismatchCount);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((int64_t)2, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(58U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((58<<4)|Cigar::ALIGN), cigarBuffer[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((42<<4)|Cigar::SOFT_CLIP), cigarBuffer[3]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)42, fragments[1][0].mismatchCount);
}

void TestFragmentBuilder::testTrailingSoftClips()
{
//    const unsigned s0 = 0;
//    const unsigned s1 = 3;
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(contigList.contigBeginOffset(4) + 16, false));
    matchList.push_back(Match(contigList.contigBeginOffset(4) + 10, true));
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 1, cluster4t, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 1, matchList.end(), cluster4t, true, fragments[1]);
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((int64_t)16, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(44U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((44<<4)|Cigar::ALIGN), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((56<<4)|Cigar::SOFT_CLIP), cigarBuffer[1]);
//    CPPUNIT_ASSERT_EQUAL((unsigned)56, fragments[0][0].mismatchCount);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].mismatchCount);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((int64_t)0, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(42U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)2, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((58<<4)|Cigar::SOFT_CLIP), cigarBuffer[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((42<<4)|Cigar::ALIGN), cigarBuffer[3]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)58, fragments[1][0].mismatchCount);
}

void TestFragmentBuilder::testLeadingAndTrailingSoftClips()
{
//    const unsigned s0 = 1;
//    const unsigned s1 = 4;
    // build fragments for a single seed on each read
    // ReferencePosition(contigId, position)
    using isaac::reference::ReferencePosition;
    using isaac::alignment::Match;
    matchList.push_back(Match(contigList.contigBeginOffset(4) + 22, false));
    matchList.push_back(Match(contigList.contigBeginOffset(4) + 11, true));
    using isaac::alignment::templateBuilder::FragmentBuilder;
    using isaac::alignment::Cigar;
    isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    isaac::alignment::Cigar cigarBuffer;
    std::vector<isaac::alignment::FragmentMetadataList> fragments(2);
    FragmentBuilder fragmentBuilder(true, flowcells, 123, 16, 1234, 3, 8, 2, false, 32, false, true, alignmentCfg, cigarBuffer, false);
    // build the fragments
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[0], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin(), matchList.begin() + 1, cluster4lt, true, fragments[0]);
//    fragmentBuilder.build(contigList, contigAnnotations, readMetadataList[1], seedMetadataList, testAdapters,
//                          isaac::alignment::TemplateLengthStatistics(), matchList.begin() + 1, matchList.end(), cluster4lt, true, fragments[1]);
    // Fragment for the first read (forward)
    CPPUNIT_ASSERT_EQUAL((int64_t)0, fragments[0][0].position);
    CPPUNIT_ASSERT_EQUAL(60U, fragments[0][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )0, fragments[0][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)false, fragments[0][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragments[0][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((10<<4)|Cigar::SOFT_CLIP), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((60<<4)|Cigar::ALIGN), cigarBuffer[1]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((30<<4)|Cigar::SOFT_CLIP), cigarBuffer[2]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[0][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)40, fragments[0][0].mismatchCount);
    // Fragment for the second read (reverse)
    CPPUNIT_ASSERT_EQUAL((int64_t)0, fragments[1][0].position);
    CPPUNIT_ASSERT_EQUAL(60U, fragments[1][0].getObservedLength());
    CPPUNIT_ASSERT_EQUAL((unsigned )1, fragments[1][0].readIndex);
    CPPUNIT_ASSERT_EQUAL((bool)true, fragments[1][0].reverse);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragments[1][0].cigarOffset);
    CPPUNIT_ASSERT_EQUAL((unsigned)3, fragments[1][0].cigarLength);
    //CPPUNIT_ASSERT_EQUAL((unsigned)((100<<4)|Cigar::MATCH), cigarBuffer[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((25<<4)|Cigar::SOFT_CLIP), cigarBuffer[3]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((60<<4)|Cigar::ALIGN), cigarBuffer[4]);
    CPPUNIT_ASSERT_EQUAL((unsigned)((15<<4)|Cigar::SOFT_CLIP), cigarBuffer[5]);
    CPPUNIT_ASSERT_EQUAL((unsigned)0, fragments[1][0].mismatchCount);
//    CPPUNIT_ASSERT_EQUAL((unsigned)40, fragments[1][0].mismatchCount);
}

