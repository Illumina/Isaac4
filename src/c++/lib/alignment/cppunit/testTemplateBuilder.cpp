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

using namespace std;

#include "RegistryName.hh"
#include "testTemplateBuilder.hh"
#include "BuilderInit.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestTemplateBuilder, registryName("TemplateBuilder"));

isaac::alignment::FragmentMetadata getFragmentMetadata(
    unsigned contigId,
    int64_t position,
    int64_t observedLength,
    unsigned readIndex,
    bool reverse,
    unsigned cigarOffset,
    unsigned cigarLength,
    const isaac::alignment::Cigar *cigarBuffer,
    unsigned mismatchCount,
    double logProbability,
    unsigned alignmentScore,
    const isaac::alignment::Cluster *cluster)
{
    isaac::alignment::FragmentMetadata f;
    f.contigId = contigId;
    f.position = position;
    f.rStrandPos = isaac::reference::ReferencePosition(contigId, position + observedLength);
    f.readIndex = readIndex;
    f.reverse = reverse;
    f.cigarOffset = cigarOffset;
    f.cigarLength = cigarLength;
    f.cigarBuffer = cigarBuffer;
    f.mismatchCount = mismatchCount;
    f.logProbability = logProbability;
    f.alignmentScore = alignmentScore;
    f.cluster = cluster;
    return f;
}

TestTemplateBuilder::TestTemplateBuilder()
    : readMetadataList(getReadMetadataList())
    , flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, isaac::flowcell::FastqFlowcellData(false, '!', false), 8, 0, std::vector<unsigned>(),
                                          readMetadataList, "blah"))
    , contigList(getContigList())
    , restOfGenomeCorrection(contigList, readMetadataList)
    , cigarBuffer(1000, 1600)
    , tls()
    , bcl0(getBcl(readMetadataList, contigList, 0, 2, 3))
    , bcl2(getBcl(readMetadataList, contigList, 2, 1, 2))
    , tile0(32)
    , tile2(31)
    , clusterId0(1234)
    , clusterId2(1234)
    , cluster0(isaac::flowcell::getMaxReadLength(readMetadataList))
    , cluster2(isaac::flowcell::getMaxReadLength(readMetadataList))
      //, f0_0(getFragmentMetadata(0,2,100,0, false, 0, 1, &cigarBuffer, 0, -78.0, 3, 254, &cluster0))
      //, f0_1(getFragmentMetadata(0,107,99,1, true, 1, 1, &cigarBuffer, 2, -92.0, 1, 253, &cluster0))
      , f0_0(getFragmentMetadata(0,2,100,0, false, 0, 1, &cigarBuffer, 0, -8.0, 254, &cluster0))
      , f0_1(getFragmentMetadata(0,107,99,1, true, 1, 1, &cigarBuffer, 2, -12.0, 253, &cluster0))
      , f1_0(getFragmentMetadata(0,7,100,0, false, 0, 1, &cigarBuffer, 0, -8.0, 254, &cluster0))
      , f1_1(getFragmentMetadata(0,0,99,1, true, 1, 1, &cigarBuffer, 2, -12.0, 253, &cluster0))
{
    cluster0.init(readMetadataList, bcl0.cluster(0), tile0, clusterId0, isaac::alignment::ClusterXy(0,0), true, 0, 0);
    cluster2.init(readMetadataList, bcl2.cluster(0), tile2, clusterId2, isaac::alignment::ClusterXy(0,0), true, 0, 0);
}

void TestTemplateBuilder::setUp()
{
}

void TestTemplateBuilder::tearDown()
{
}

void TestTemplateBuilder::testConstructor()
{
    CPPUNIT_ASSERT_EQUAL(150U, tls.getMin());
    CPPUNIT_ASSERT_EQUAL(190U, tls.getMedian());
    CPPUNIT_ASSERT_EQUAL(250U, tls.getMax());
    CPPUNIT_ASSERT_EQUAL(20U, tls.getLowStdDev());
    CPPUNIT_ASSERT_EQUAL(30U, tls.getHighStdDev());
    CPPUNIT_ASSERT(tls.isStable());
    CPPUNIT_ASSERT_EQUAL(std::string("FR+"), tls.alignmentModelName(tls.getBestModel(0)));
    CPPUNIT_ASSERT_EQUAL(std::string("RF-"), tls.alignmentModelName(tls.getBestModel(1)));
    CPPUNIT_ASSERT_EQUAL(1U, readMetadataList[1].getIndex());
    CPPUNIT_ASSERT_EQUAL(1U, cluster0[1].getIndex());
    CPPUNIT_ASSERT_EQUAL(254U, f0_0.alignmentScore);
}

void TestTemplateBuilder::checkUnalignedTemplate(
    const isaac::alignment::BamTemplate &bamTemplate,
    const isaac::alignment::Cluster &cluster) const
{
    checkUnalignedFragment(bamTemplate, cluster, 0, 0);
    checkUnalignedFragment(bamTemplate, cluster, 1, 1);
}


void TestTemplateBuilder::checkUnalignedFragment(
    const isaac::alignment::BamTemplate &bamTemplate,
    const isaac::alignment::Cluster &cluster,
    const unsigned i,
    const unsigned readIndex) const
{
    CPPUNIT_ASSERT_EQUAL(i, readIndex);
    // check for an unaligned fragment
    CPPUNIT_ASSERT(bamTemplate.getFragmentMetadata(i).isNoMatch());
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).getObservedLength());
    CPPUNIT_ASSERT_EQUAL(readIndex, bamTemplate.getFragmentMetadata(i).readIndex);
    CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(i).reverse);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).cigarLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(i).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(0.0, bamTemplate.getFragmentMetadata(i).logProbability);
    CPPUNIT_ASSERT_EQUAL(-1U, bamTemplate.getFragmentMetadata(i).alignmentScore);
    CPPUNIT_ASSERT_EQUAL(&cluster, bamTemplate.getFragmentMetadata(i).cluster);
}

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;

void TestTemplateBuilder::testEmptyMatchList()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    const isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    std::auto_ptr<TemplateBuilder> templateBuilder(new TemplateBuilder(true, flowcells, 10, 10, 16, 4, 1000, 1000, 1000, false, true, false, false, 8, 2, false, 32, false,
                                                                       alignmentCfg,
                                                                       TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED, 4, false));
    BamTemplate bamTemplate;
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentCount());
    isaac::alignment::TemplateBuilder::FragmentMetadataLists fragments;
    templateBuilder->buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentCount());
    checkUnalignedTemplate(bamTemplate, cluster0);
    CPPUNIT_ASSERT_EQUAL(-1U, bamTemplate.getAlignmentScore());
    // initialize the first fragment with garbage
    fragments[0].push_back(f0_0);
    fragments[1].push_back(f0_1);
    templateBuilder->buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
    // clear the fragments and check again
    fragments[0].clear();
    fragments[1].clear();
    templateBuilder->buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentCount());
    checkUnalignedTemplate(bamTemplate, cluster0);
    CPPUNIT_ASSERT_EQUAL(-1U, bamTemplate.getAlignmentScore());
}

void TestTemplateBuilder::testOrphan()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    const isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    TemplateBuilder templateBuilder(true, flowcells, 10, 10, 16, 4, 1000, 1000, 1000, false, true, false, false, 8, 2, false, 32, true,
                                    alignmentCfg,
                                    TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED, 4, false);
    BamTemplate bamTemplate;
    isaac::alignment::TemplateBuilder::FragmentMetadataLists fragments;
    // align on the first read only
    fragments[0].push_back(f0_0);
    templateBuilder.buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
    // this orphan should be rescued
//    CPPUNIT_ASSERT_EQUAL(1136U, bamTemplate.getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).contigId);
    CPPUNIT_ASSERT_EQUAL(2L, bamTemplate.getFragmentMetadata(0).position);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).readIndex);
    CPPUNIT_ASSERT_EQUAL(100U, bamTemplate.getFragmentMetadata(0).getObservedLength());
    CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(0).reverse);
//    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(0).cigarLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_0.logProbability, bamTemplate.getFragmentMetadata(0).logProbability);
    CPPUNIT_ASSERT_EQUAL(534U, bamTemplate.getFragmentMetadata(0).alignmentScore);
//    CPPUNIT_ASSERT_EQUAL(569U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
    // align on the second read only
    fragments[0].clear();
    fragments[1].push_back(f0_1);
    templateBuilder.buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
    // this one should be rescued as well
//    CPPUNIT_ASSERT_EQUAL(1119U, bamTemplate.getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(1).contigId);
    CPPUNIT_ASSERT_EQUAL(107L, bamTemplate.getFragmentMetadata(1).position);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).readIndex);
    CPPUNIT_ASSERT_EQUAL(99U, bamTemplate.getFragmentMetadata(1).getObservedLength());
    CPPUNIT_ASSERT_EQUAL(true, bamTemplate.getFragmentMetadata(1).reverse);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarLength);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(1).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_1.logProbability, bamTemplate.getFragmentMetadata(1).logProbability);
    CPPUNIT_ASSERT_EQUAL(517U, bamTemplate.getFragmentMetadata(1).alignmentScore);
//    CPPUNIT_ASSERT_EQUAL(569U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);
}

void TestTemplateBuilder::testUnique()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    const isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    TemplateBuilder templateBuilder(true, flowcells, 10, 10, 16, 4, 1000, 1000, 1000, false, true, false, false, 8, 2, false, 32, true,
                                    alignmentCfg,
                                    TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED, 4, false);
    BamTemplate bamTemplate;
    isaac::alignment::TemplateBuilder::FragmentMetadataLists fragments;
    fragments[0].push_back(f0_0);
    fragments[1].push_back(f0_1);
    templateBuilder.buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
    CPPUNIT_ASSERT_EQUAL(1084U, bamTemplate.getAlignmentScore());
    //CPPUNIT_ASSERT_EQUAL(bamTemplate.getFragmentMetadata(0).getAlignmentScore() + bamTemplate.getFragmentMetadata(1).getAlignmentScore(),
    //                     bamTemplate.getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(534U, bamTemplate.getFragmentMetadata(0).getAlignmentScore());
    CPPUNIT_ASSERT_EQUAL(517U, bamTemplate.getFragmentMetadata(1).getAlignmentScore());
    // check the first read
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).contigId);
    CPPUNIT_ASSERT_EQUAL(2L, bamTemplate.getFragmentMetadata(0).position);
    CPPUNIT_ASSERT_EQUAL(100U, bamTemplate.getFragmentMetadata(0).getObservedLength());
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).readIndex);
    CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(0).reverse);
//    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(0).cigarLength);
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_0.logProbability, bamTemplate.getFragmentMetadata(0).logProbability);
    CPPUNIT_ASSERT_EQUAL(534U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
    // check the second read
    CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(1).contigId);
    CPPUNIT_ASSERT_EQUAL(107L, bamTemplate.getFragmentMetadata(1).position);
    CPPUNIT_ASSERT_EQUAL(99U, bamTemplate.getFragmentMetadata(1).getObservedLength());
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).readIndex);
    CPPUNIT_ASSERT_EQUAL(true, bamTemplate.getFragmentMetadata(1).reverse);
//    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarOffset);
    CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarLength);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(1).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(f0_1.logProbability, bamTemplate.getFragmentMetadata(1).logProbability);
    CPPUNIT_ASSERT_EQUAL(517U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);
}

void TestTemplateBuilder::testPeAdapterTrim()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    const isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);

    // not trimming
    {
        TemplateBuilder templateBuilder(true, flowcells, 10, 10, 16, 4, 1000, 1000, 1000, false, true, false, false, 8, 2, false, 32, true,
                                        alignmentCfg,
                                        TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED, 4, false);
        BamTemplate bamTemplate;
        isaac::alignment::TemplateBuilder::FragmentMetadataLists fragments;
        fragments[0].push_back(f1_0);
        fragments[0].back().updateAlignment(
            false, alignmentCfg, readMetadataList[0], contigList,
            fragments[0].back().reverse, fragments[0].back().getContigId(), fragments[0].back().getPosition(),
            *fragments[0].back().cigarBuffer, fragments[0].back().cigarOffset, fragments[0].back().cigarLength);

        fragments[1].push_back(f1_1);
        fragments[1].back().updateAlignment(
            false, alignmentCfg, readMetadataList[1], contigList,
            fragments[1].back().reverse, fragments[1].back().getContigId(), fragments[1].back().getPosition(),
            *fragments[1].back().cigarBuffer, fragments[1].back().cigarOffset, fragments[1].back().cigarLength);

        templateBuilder.buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
        // check the first read
        CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).contigId);
        CPPUNIT_ASSERT_EQUAL(7L, bamTemplate.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(100U, bamTemplate.getFragmentMetadata(0).getObservedLength());
        CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).readIndex);
        CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(0).reverse);
        CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(0).cigarLength);
        CPPUNIT_ASSERT_EQUAL(70U, bamTemplate.getFragmentMetadata(0).mismatchCount);
        CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
        // check the second read
        CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(1).contigId);
        CPPUNIT_ASSERT_EQUAL(0L, bamTemplate.getFragmentMetadata(1).position);
        CPPUNIT_ASSERT_EQUAL(100U, bamTemplate.getFragmentMetadata(1).getObservedLength());
        CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).readIndex);
        CPPUNIT_ASSERT_EQUAL(true, bamTemplate.getFragmentMetadata(1).reverse);
        CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).cigarLength);
        CPPUNIT_ASSERT_EQUAL(86U, bamTemplate.getFragmentMetadata(1).mismatchCount);
        CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);
    }

    // trimming
    {
        TemplateBuilder templateBuilder(true, flowcells, 10, 10, 16, 4, 1000, 1000, 1000, false, true, true, false, 8, 2, false, 32, true,
                                        alignmentCfg,
                                        TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED, 4, false);
        BamTemplate bamTemplate;
        isaac::alignment::TemplateBuilder::FragmentMetadataLists fragments;
        fragments[0].push_back(f1_0);
        fragments[0].back().updateAlignment(
            false, alignmentCfg, readMetadataList[0], contigList,
            fragments[0].back().reverse, fragments[0].back().getContigId(), fragments[0].back().getPosition(),
            *fragments[0].back().cigarBuffer, fragments[0].back().cigarOffset, fragments[0].back().cigarLength);

        fragments[1].push_back(f1_1);
        fragments[1].back().updateAlignment(
            false, alignmentCfg, readMetadataList[1], contigList,
            fragments[1].back().reverse, fragments[1].back().getContigId(), fragments[1].back().getPosition(),
            *fragments[1].back().cigarBuffer, fragments[1].back().cigarOffset, fragments[1].back().cigarLength);

        templateBuilder.buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);
        // check the first read
        CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).contigId);
        CPPUNIT_ASSERT_EQUAL(7L, bamTemplate.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(93U, bamTemplate.getFragmentMetadata(0).getObservedLength());
        CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(0).readIndex);
        CPPUNIT_ASSERT_EQUAL(false, bamTemplate.getFragmentMetadata(0).reverse);
        CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(0).cigarLength);
        CPPUNIT_ASSERT_EQUAL(67U, bamTemplate.getFragmentMetadata(0).mismatchCount);
        CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
        // check the second read
        CPPUNIT_ASSERT_EQUAL(0U, bamTemplate.getFragmentMetadata(1).contigId);
        CPPUNIT_ASSERT_EQUAL(7L, bamTemplate.getFragmentMetadata(1).position);
        CPPUNIT_ASSERT_EQUAL(93U, bamTemplate.getFragmentMetadata(1).getObservedLength());
        CPPUNIT_ASSERT_EQUAL(1U, bamTemplate.getFragmentMetadata(1).readIndex);
        CPPUNIT_ASSERT_EQUAL(true, bamTemplate.getFragmentMetadata(1).reverse);
        CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(1).cigarLength);
        CPPUNIT_ASSERT_EQUAL(79U, bamTemplate.getFragmentMetadata(1).mismatchCount);
        CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);
    }
}

/**
 * \brief This test was originally designed to ensure the pair that matches the tls is picked.
 * Currently, everything on the same contig with the correct orientation with size below max_ + 50000
 * is considered allright. So, verify that the best alignment score pair is picked.
 */
void TestTemplateBuilder::testMultiple()
{
    using isaac::alignment::TemplateBuilder;
    using isaac::alignment::BamTemplate;
    using isaac::alignment::FragmentMetadata;
    using isaac::alignment::BandedSmithWaterman;
    const isaac::alignment::AlignmentCfg alignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, 20000);
    TemplateBuilder templateBuilder(true, flowcells, 10, 10, 16, 4, 1000, 1000, 1000, false, true, false, false, 8, 2, false, 32, true,
                                    alignmentCfg,
                                    TemplateBuilder::DODGY_ALIGNMENT_SCORE_UNALIGNED, 4, false);
    BamTemplate bamTemplate;

    isaac::alignment::TemplateBuilder::FragmentMetadataLists fragments;
    FragmentMetadata t0 = f0_0;
    FragmentMetadata t1 = f0_1;
    for (size_t i = 0; 2 > i; ++i)
    {
        fragments[0].push_back(t0);
        t0.position += 56;
        t0.rStrandPos += 56;
        fragments[0].push_back(t0);
        t0.position += 65;
        t0.rStrandPos += 65;
        fragments[1].push_back(t1);
        t1.position += 300;
        t1.rStrandPos += 300;
    }
    t0 = f0_0;
    t1 = f0_1;
    t0.contigId = 1;
    t0.rStrandPos = isaac::reference::ReferencePosition(t0.contigId, t0.rStrandPos.getPosition());
    t1.contigId = 1;
    t1.rStrandPos = isaac::reference::ReferencePosition(t1.contigId, t1.rStrandPos.getPosition());
    for (size_t i = 0; 2 > i; ++i)
    {
        t0.position += 56;
        t0.rStrandPos += 56;
        fragments[0].push_back(t0);
        t0.position += 65;
        t0.rStrandPos += 65;
        fragments[0].push_back(t0);
        t1.position += 401;
        t1.rStrandPos += 401;
        fragments[1].push_back(t1);
    }
    t0 = f0_0;
    t1 = f0_1;
    t0.contigId = 1;
    t0.rStrandPos = isaac::reference::ReferencePosition(t0.contigId, t0.rStrandPos.getPosition());
    t1.contigId = 1;
    t1.rStrandPos = isaac::reference::ReferencePosition(t1.contigId, t1.rStrandPos.getPosition());
    t0.logProbability += 2;
    t1.logProbability += 2;
    fragments[0].push_back(t0);
    fragments[1].push_back(t1);
    t0.logProbability -= 2;
    t1.logProbability -= 2;
    for (size_t i = 0; 2 > i; ++i)
    {
        t0.position += 36;
        t0.rStrandPos += 36;
        fragments[0].push_back(t0);
        t0.position += 45;
        t0.rStrandPos += 45;
        fragments[0].push_back(t0);
        t1.position += 402;
        t1.rStrandPos += 402;
        fragments[1].push_back(t1);
    }
    std::sort(fragments[0].begin(), fragments[0].end(), FragmentMetadata::bestUngappedLess);
    std::sort(fragments[1].begin(), fragments[1].end(), FragmentMetadata::bestUngappedLess);
    FragmentMetadata best0 = fragments[0].front();
    FragmentMetadata best1 = fragments[1].front();
    templateBuilder.buildCombinationTemplate(contigList, restOfGenomeCorrection, readMetadataList, fragments, cluster0, tls, bamTemplate);

    CPPUNIT_ASSERT_EQUAL(1101U, bamTemplate.getAlignmentScore());
    // check the first read
    CPPUNIT_ASSERT_EQUAL(unsigned(best0.getFStrandReferencePosition().getContigId()), bamTemplate.getFragmentMetadata(0).contigId);
    CPPUNIT_ASSERT_EQUAL(int64_t(best0.getFStrandReferencePosition().getPosition()), bamTemplate.getFragmentMetadata(0).position);
    CPPUNIT_ASSERT_EQUAL(best0.getObservedLength(), unsigned(bamTemplate.getFragmentMetadata(0).getObservedLength()));
    CPPUNIT_ASSERT_EQUAL(best0.getReadIndex(), bamTemplate.getFragmentMetadata(0).readIndex);
    CPPUNIT_ASSERT_EQUAL(best0.isReverse(), bamTemplate.getFragmentMetadata(0).reverse);
    CPPUNIT_ASSERT_EQUAL(best0.getCigarLength(), bamTemplate.getFragmentMetadata(0).cigarLength);
    CPPUNIT_ASSERT_EQUAL(best0.getMismatchCount(), bamTemplate.getFragmentMetadata(0).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(best0.logProbability, bamTemplate.getFragmentMetadata(0).logProbability);
    CPPUNIT_ASSERT_EQUAL(2U, bamTemplate.getFragmentMetadata(0).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(0).cluster);
    // check the second read
    CPPUNIT_ASSERT_EQUAL(unsigned(best1.getFStrandReferencePosition().getContigId()), bamTemplate.getFragmentMetadata(1).contigId);
    CPPUNIT_ASSERT_EQUAL(int64_t(best1.getFStrandReferencePosition().getPosition()), bamTemplate.getFragmentMetadata(1).position);
    CPPUNIT_ASSERT_EQUAL(best1.getObservedLength(), unsigned(bamTemplate.getFragmentMetadata(1).getObservedLength()));
    CPPUNIT_ASSERT_EQUAL(best1.getReadIndex(), bamTemplate.getFragmentMetadata(1).readIndex);
    CPPUNIT_ASSERT_EQUAL(best1.isReverse(), bamTemplate.getFragmentMetadata(1).reverse);
    CPPUNIT_ASSERT_EQUAL(best1.getCigarLength(), bamTemplate.getFragmentMetadata(1).cigarLength);
    CPPUNIT_ASSERT_EQUAL(best1.getMismatchCount(), bamTemplate.getFragmentMetadata(1).mismatchCount);
    CPPUNIT_ASSERT_EQUAL(best1.logProbability, bamTemplate.getFragmentMetadata(1).logProbability);
    CPPUNIT_ASSERT_EQUAL(3U, bamTemplate.getFragmentMetadata(1).alignmentScore);
    CPPUNIT_ASSERT(&cluster0 == bamTemplate.getFragmentMetadata(1).cluster);

}

void TestTemplateBuilder::testAll()
{
    {
        testConstructor();
    //    testEmptyMatchList();
        testOrphan();
        testUnique();
        testPeAdapterTrim();
        testMultiple();
    }
}

DummyTemplateLengthStatistics::DummyTemplateLengthStatistics():
    TemplateLengthStatistics()
{
    setMin(150, -1);
    setMax(250, -1);
    setMedian(190, -1);
    setLowStdDev(20);
    setHighStdDev(30);
    setBestModel(FRp, 0); // FR+
    setBestModel(RFm, 1); // RF-
    setStable(true);
}
