/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file testFragmentBuilder2.cpp
 **
 ** More fragment builder tests.
 **
 ** \author Roman Petrovski
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/lambda/construct.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testFragmentBuilder2.hh"

#include "alignment/templateBuilder/GappedAligner.hh"
#include "alignment/templateBuilder/UngappedAligner.hh"
#include "alignment/templateBuilder/FragmentSequencingAdapterClipper.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cluster.hh"
#include "flowcell/SequencingAdapterMetadata.hh"

#include "BuilderInit.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestFragmentBuilder2, registryName("FragmentBuilder2"));

static const isaac::alignment::SequencingAdapterList noAdapters;

TestFragmentBuilder2::TestFragmentBuilder2() :
    readMetadataList(getReadMetadataList()),
    flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, isaac::flowcell::FastqFlowcellData(false, '!', false), 8, 0, std::vector<unsigned>(),
                                         readMetadataList, "blah"))
{
    cigarBuffer_.reserve(1024);
}

void TestFragmentBuilder2::setUp()
{
}

void TestFragmentBuilder2::tearDown()
{
}

namespace testFragmentBuilder2
{

struct ReadInit : public std::pair<std::string, std::string>
{
    static std::string reverseString(std::string fwd)
    {
        std::reverse(fwd.begin(), fwd.end());
        return fwd;
    }

    typedef std::pair<std::string, std::string> BaseType;
    ReadInit(const std::string &read, const std::string &qual, const bool reverse) :
        BaseType(reverse ? reverseString(read) : read, qual)
    {

    }
};

} // namespace testFragmentBuilder2

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> testFragmentBuilder2::ReadInit& operator >><testFragmentBuilder2::ReadInit >(
    testFragmentBuilder2::ReadInit &input,
    isaac::alignment::Read &read)
{
    ISAAC_ASSERT_MSG(input.first.length() == input.second.length(), "sequence and quality must be of equal lengths");

    read.forwardSequence_ = vectorFromString(input.first);
    read.forwardQuality_ = vectorFromString(input.second);

    std::for_each(read.forwardQuality_.begin(), read.forwardQuality_.end(), &phredToBcl);

    read.reverseSequence_ = read.forwardSequence_;
    read.reverseQuality_ = read.forwardQuality_;
    std::reverse(read.reverseSequence_.begin(), read.reverseSequence_.end());
    std::reverse(read.reverseQuality_.begin(), read.reverseQuality_.end());

    return input;
}

}
}

void TestFragmentBuilder2::testEverything()
{
    ISAAC_SCOPE_BLOCK_CERR
    {
    testMismatchCount();
    testMismatchCycles();
    testMismatchCyclesWithSoftClip();
    testGapped();
    testGappedWithNs();
    }
}

static const int ELAND_MATCH_SCORE = 2;
static const int ELAND_MISMATCH_SCORE = -1;
static const int ELAND_GAP_OPEN_SCORE = -15;
static const int ELAND_GAP_EXTEND_SCORE = -3;
static const int ELAND_MIN_GAP_EXTEND_SCORE = 25;


void TestFragmentBuilder2::align(
    const std::string &read,
    const std::string &reference,
    const isaac::alignment::SequencingAdapterList &adapters,
    isaac::alignment::FragmentMetadata &fragmentMetadata,
    const bool gapped)
{
    static const std::string irrelevantQualities("CFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGGCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGGCFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGG");

    align(read, irrelevantQualities.substr(0, read.length()), reference, adapters, fragmentMetadata, gapped);
}

void TestFragmentBuilder2::align(
    const std::string &read,
    const std::string &qual,
    const std::string &reference,
    const isaac::alignment::SequencingAdapterList &adapters,
    isaac::alignment::FragmentMetadata &fragmentMetadata,
    const bool gapped)
{
    align(read, qual, reference, adapters, fragmentMetadata,
          isaac::alignment::AlignmentCfg(ELAND_MATCH_SCORE, ELAND_MISMATCH_SCORE, ELAND_GAP_OPEN_SCORE, ELAND_GAP_EXTEND_SCORE, ELAND_MIN_GAP_EXTEND_SCORE, -1U), gapped);
}

void TestFragmentBuilder2::align(
    const std::string &read,
    const std::string &qual,
    const std::string &reference,
    const isaac::alignment::SequencingAdapterList &adapters,
    isaac::alignment::FragmentMetadata &fragmentMetadata,
    const isaac::alignment::AlignmentCfg &alignmentCfg,
    const bool gapped,
    const unsigned int endCyclesMasked)
{
    isaac::alignment::Cluster cluster(isaac::flowcell::getMaxReadLength(flowcells));
    testFragmentBuilder2::ReadInit init(read, qual, fragmentMetadata.reverse);
    init >> cluster.at(0);
    cluster.at(0).maskCyclesFromEnd(endCyclesMasked);

    if (fragmentMetadata.isNoMatch())
    {
        fragmentMetadata.contigId = 0;
        fragmentMetadata.position = 0;
    }
    fragmentMetadata.cluster = &cluster;
    fragmentMetadata.cigarBuffer = &cigarBuffer_;

    TestContigList contigList(reference);
    const isaac::reference::Contig &referenceContig = contigList.at(0);

    isaac::alignment::templateBuilder::FragmentSequencingAdapterClipper adapterClipper(adapters);
    adapterClipper.checkInitStrand(fragmentMetadata, referenceContig);

    isaac::alignment::templateBuilder::UngappedAligner ungappedAligner(true, alignmentCfg);

    ungappedAligner.alignUngapped(fragmentMetadata, cigarBuffer_, readMetadataList[fragmentMetadata.getReadIndex()], adapterClipper, contigList);
    if (gapped)
    {
        isaac::alignment::templateBuilder::GappedAligner gappedAligner(true, flowcells, false, 32, alignmentCfg);
        isaac::alignment::FragmentMetadata tmp = fragmentMetadata;
        const unsigned matchCount = gappedAligner.alignGapped(
            readMetadataList[fragmentMetadata.getReadIndex()], adapterClipper, contigList, tmp, cigarBuffer_);
        if (matchCount + isaac::alignment::BandedSmithWaterman<16>::WIDEST_GAP_SIZE > fragmentMetadata.getObservedLength() &&
                                (tmp.mismatchCount <= 5) &&
                                (fragmentMetadata.mismatchCount > tmp.mismatchCount) &&
                                fragmentMetadata.logProbability < tmp.logProbability)
        {
            fragmentMetadata = tmp;
        }
    }

    fragmentMetadata.recomputeHeadAnchor(contigList);
    fragmentMetadata.recomputeTailAnchor(contigList);
}


void TestFragmentBuilder2::testMismatchCount()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    align("TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //||||||||||||||||||||||||||||||||||||||||x|||xx|x|xx|xxxxx|x|xxxxx|x||x|x|x|x||||||||||x|x|x|x|xxx|x|
          "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTCTCTTCTCTGGAATATGATAAAAAAAAAAAAAAAAAGTGCACCGCCAAAAAAAAAAAAAA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("100M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(30U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(30U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(100U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 99U), fragmentMetadata.getRStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(100,100, false), fragmentMetadata.tailAnchor());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(0,16, false), fragmentMetadata.headAnchor());
}

void TestFragmentBuilder2::testMismatchCycles()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = true;

    align("TGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         //||||||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("100M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(100U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getFStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(92U, *fragmentMetadata.getMismatchCyclesBegin());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(0,16, false), fragmentMetadata.tailAnchor());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(84,100, false), fragmentMetadata.headAnchor());
}

void TestFragmentBuilder2::testMismatchCyclesWithSoftClip()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    fragmentMetadata.contigId = 0;
    fragmentMetadata.position = -2;

    align("TTTGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTA",
         //**||||||||x|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            "TGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata);

    CPPUNIT_ASSERT_EQUAL(std::string("2S98M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(98U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 0U), fragmentMetadata.getFStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(11U, *fragmentMetadata.getMismatchCyclesBegin());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(84,100, false), fragmentMetadata.tailAnchor());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(2,18, false), fragmentMetadata.headAnchor());
}

void TestFragmentBuilder2::testGapped()
{
    isaac::alignment::FragmentMetadata fragmentMetadata;
    fragmentMetadata.reverse = false;
    fragmentMetadata.contigId = 0;
    fragmentMetadata.position = 1;

    align( "TTTGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTA",
         // ||||||||||x||||||||||||||||||||||||||||||||||||||||||^\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
          "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
         noAdapters,
         fragmentMetadata, true);

    CPPUNIT_ASSERT_EQUAL(std::string("53M1D47M"), fragmentMetadata.getCigarString());
    CPPUNIT_ASSERT_EQUAL(1U, fragmentMetadata.getMismatchCount());
    CPPUNIT_ASSERT_EQUAL(2U, fragmentMetadata.getEditDistance());
    CPPUNIT_ASSERT_EQUAL(101U, fragmentMetadata.getObservedLength());
    CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 1U), fragmentMetadata.getFStrandReferencePosition());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(84,100, false), fragmentMetadata.tailAnchor());
    CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(0,16, false), fragmentMetadata.headAnchor());
}

void TestFragmentBuilder2::testGappedWithNs()
{
    {
        isaac::alignment::FragmentMetadata fragmentMetadata;
        fragmentMetadata.reverse = false;
        fragmentMetadata.contigId = 0;
        fragmentMetadata.position = 1;

        align( "TTTGGTTAAGATAGCGGTAAAAGCGTGTTACCGCAATGTTCTGNNNNTTATACACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTA",
             // |||||||||||||||||||||||||||||||||||||||||||||||||||||^\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
             // ||||||||x|||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||
              "ATTTGGTTAAGGTAGCGGTAAAAGCGTGTTACCGCAATGTTCTGTCTCTTATACAACATCTAGATGTGTATAAGAGACAGGTGCACCGCCTATACACATCTAGAATAAGAGACAGGTGCACCGCCTATACACATCTAGA",
             noAdapters,
             fragmentMetadata, true);

        CPPUNIT_ASSERT_EQUAL(std::string("53M1D47M"), fragmentMetadata.getCigarString());
        CPPUNIT_ASSERT_EQUAL(5U, fragmentMetadata.getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(6U, fragmentMetadata.getEditDistance());
        CPPUNIT_ASSERT_EQUAL(101U, fragmentMetadata.getObservedLength());
        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 1U), fragmentMetadata.getFStrandReferencePosition());
        CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(84,100, false), fragmentMetadata.tailAnchor());
        CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(0,16, false), fragmentMetadata.headAnchor());
    }

    {
        isaac::alignment::FragmentMetadata fragmentMetadata;
        fragmentMetadata.reverse = true;
        fragmentMetadata.contigId = 0;
        fragmentMetadata.position = 1;

        align( "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAATGGAATCAAAATAAAAAGGNNNCAAACGGAATTATCGAA",
               "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DGGGGFIGCGC99E?<C<++33!!!HJIHHHDHFADDAB1@?",
              "AAAACGGAATCAAATGGAATTATCAAATGCAATCGAAGAGAATCATCGAATGATGGACTCAAATGGAATCAACGTCAAACGGAATCAAATGGAATTATCAAATGCAATCGAAGAGAATCATCGAATGGACTCGAATGGAACCATCTAATGGAATGGAATGGAATAATCCATGGACTCGAATGCAATCATCATCAAATGGAATCGAATGGAATCATCGAATGGACTCAAATGGAATAATCATTGAACGGAATCAAATGGAATCATCATCGGATGGAA",
             noAdapters,
             fragmentMetadata, isaac::alignment::AlignmentCfg(0, -4, -6, -1, 20, -1U), true, 58);

        CPPUNIT_ASSERT_EQUAL(std::string("58S43M"), fragmentMetadata.getCigarString());
        CPPUNIT_ASSERT_EQUAL(11U, fragmentMetadata.getMismatchCount());
        CPPUNIT_ASSERT_EQUAL(11U, fragmentMetadata.getEditDistance());
        CPPUNIT_ASSERT_EQUAL(43U, fragmentMetadata.getObservedLength());
        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 59U), fragmentMetadata.getFStrandReferencePosition());
        CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(58,58, false), fragmentMetadata.tailAnchor());
        CPPUNIT_ASSERT_EQUAL(isaac::alignment::Anchor(85,101, false), fragmentMetadata.headAnchor());
    }
}
