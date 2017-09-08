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
 ** \file testOverlappingEndsClipper.cpp
 **
 ** \author Roman Petrovski
 **/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/assign.hpp>

#include "alignment/Cluster.hh"
#include "alignment/matchSelector/OverlappingEndsClipper.hh"

#include "BuilderInit.hh"
#include "RegistryName.hh"
#include "testOverlappingEndsClipper.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestOverlappingEndsClipper, registryName("OverlappingEndsClipper"));

TestOverlappingEndsClipper::TestOverlappingEndsClipper() : cluster_(1234)
{
    cigarBuffer_.reserve(1024);
}

void TestOverlappingEndsClipper::setUp()
{
}

void TestOverlappingEndsClipper::tearDown()
{
}

namespace testOverlappingEndsClipper
{

//static const std::string irrelevantQualities("CFCEEBFHEHDGBDBEDDEGEHHFHEGBHHDDDB<F>FGGBFGGFGCGGGDGGDDFHHHFEGGBGDGGBGGBEGEGGBGEHDHHHGGGGGDGGGG?GGGG");

struct ReadInit : public std::pair<std::string, std::string>
{
    static std::string reverseString(std::string fwd)
    {
        std::reverse(fwd.begin(), fwd.end());
        return fwd;
    }

    typedef std::pair<std::string, std::string> BaseType;
    ReadInit(const std::string &read, const std::string &quality, const bool reverse) :
        BaseType(reverse ? reverseString(read) : read, reverse ? quality : reverseString(quality))
    {

    }
};

} // namespace testSemialignedClipper

namespace isaac
{
namespace alignment
{

inline void phredToBcl(char &qual)
{
    qual -= 33;
}

template<class InpuT> InpuT& operator >>(InpuT &input, isaac::alignment::Read &read);

template<> testOverlappingEndsClipper::ReadInit& operator >><testOverlappingEndsClipper::ReadInit >(
    testOverlappingEndsClipper::ReadInit &input,
    isaac::alignment::Read &read)
{
//    ISAAC_THREAD_CERR << input.first << " " << input.second << std::endl;
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

void TestOverlappingEndsClipper::testEverything()
{
    ISAAC_SCOPE_BLOCK_CERR
    {
    TestContigList contigList;
    isaac::alignment::BamTemplate templ;


    {
        init("ACGT", "CFCE", false,
             "   ACGT", "BDBE", true,
             0,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(3L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("3M1S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(3L, templ.getFragmentMetadata(1).position);
    }

    {
        init("ACGT", "CFCE", false,
             " ACGT", "BDBE", true,
             0,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("3M1S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("2S2M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(3L, templ.getFragmentMetadata(1).position);
    }

    {
        init("ACGT", "BAAA", false,
             " ACGT", "CFCE", true,
             0,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("2M2S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("1S3M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(2L, templ.getFragmentMetadata(1).position);
    }

    {
        init("ACGT", "BAAA", false,
             "ACGT", "CFCE", true,
             0,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("1M3S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("1S3M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);
    }

    {// only PE pairs with adapters trimmed can be clipped. Mate-pair does not get trimmed
        init(" ACGT", "BAAA", false,
             "ACGT", "CFCE", true,
             0,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(1).position);
    }

    {// only PE pairs with adapters trimmed can be clipped. Mate-pair does not get trimmed
        init(" ACGT", "BAAA", false,
             "ACGT", "CFCE", true,
             1,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("3M1S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("1S3M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("1M3S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(1L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("2S2M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(2L, templ.getFragmentMetadata(1).position);
    }

    {// check overlap of one base special case
        init("ACGT", "BAAA", false,
             "   ACGT", "CFCE", true,
             0,
             "ACGT", templ, contigList);

        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(3L, templ.getFragmentMetadata(1).position);

        isaac::alignment::matchSelector::OverlappingEndsClipper clipper;
        clipper.clip(contigList, templ);

        CPPUNIT_ASSERT_EQUAL(std::string("3M1S"), templ.getFragmentMetadata(0).getCigarString());
        CPPUNIT_ASSERT_EQUAL(0L, templ.getFragmentMetadata(0).position);
        CPPUNIT_ASSERT_EQUAL(std::string("4M"), templ.getFragmentMetadata(1).getCigarString());
        CPPUNIT_ASSERT_EQUAL(3L, templ.getFragmentMetadata(1).position);
    }
    }
}

static TestContigList makeContigList(const std::string forward, int64_t &firstPosOffset)
{
    std::string::const_iterator begin = std::find_if(forward.begin(), forward.end(),
                                                     [](char c){return c != ' ';});

    firstPosOffset = -std::distance(forward.begin(), begin);
    return TestContigList (std::string(begin, forward.end()));
}

void TestOverlappingEndsClipper::init(
    const std::string &read1, const std::string &quality1, const bool read1Reverse,
    const std::string &read2, const std::string &quality2, const bool read2Reverse,
    const unsigned overhang,
    const std::string &reference,
    isaac::alignment::BamTemplate &templ,
    TestContigList &contigList)
{
    int64_t firstPosOffset = 0;
    contigList = makeContigList(reference, firstPosOffset);

    const unsigned r1Start = read1.find_first_not_of(' ');
    const unsigned r2Start = read2.find_first_not_of(' ');
    static isaac::flowcell::ReadMetadataList readMetadataList;
    readMetadataList = getReadMetadataList(read1.length() - r1Start, read2.length() - r2Start);

    testOverlappingEndsClipper::ReadInit init1(read1.substr(r1Start), quality1, read1Reverse);
    init1 >> cluster_.at(0);
    testOverlappingEndsClipper::ReadInit init2(read2.substr(r2Start), quality2, read2Reverse);
    init2 >> cluster_.at(1);

    templ.reset(readMetadataList, cluster_);

    cigarBuffer_.clear();

    templ.getFragmentMetadata(0).reverse = read1Reverse;
    templ.getFragmentMetadata(0).cigarBuffer = &cigarBuffer_;
    templ.getFragmentMetadata(0).cigarOffset = cigarBuffer_.size();
    templ.getFragmentMetadata(0).contigId = 0;
    templ.getFragmentMetadata(0).position = read1.find_first_not_of(' ');
    templ.getFragmentMetadata(0).rStrandPos =
        isaac::reference::ReferencePosition(0 , templ.getFragmentMetadata(0).position + readMetadataList.at(0).getLength() - overhang);
    cigarBuffer_.addOperation(readMetadataList.at(0).getLength() - overhang, isaac::alignment::Cigar::ALIGN);
    if (overhang)
    {
        cigarBuffer_.addOperation(overhang, isaac::alignment::Cigar::SOFT_CLIP);
    }
    templ.getFragmentMetadata(0).cigarLength = cigarBuffer_.size() - templ.getFragmentMetadata(0).cigarOffset;

    templ.getFragmentMetadata(1).reverse = read2Reverse;
    templ.getFragmentMetadata(1).cigarBuffer = &cigarBuffer_;
    templ.getFragmentMetadata(1).cigarOffset = cigarBuffer_.size();
    templ.getFragmentMetadata(1).contigId = 0;
    templ.getFragmentMetadata(1).position = read2.find_first_not_of(' ') + overhang;
    templ.getFragmentMetadata(1).rStrandPos =
        isaac::reference::ReferencePosition(0 , templ.getFragmentMetadata(1).position + readMetadataList.at(1).getLength() - overhang);
    if (overhang)
    {
        cigarBuffer_.addOperation(overhang, isaac::alignment::Cigar::SOFT_CLIP);
    }
    cigarBuffer_.addOperation(readMetadataList.at(1).getLength() - overhang, isaac::alignment::Cigar::ALIGN);
    templ.getFragmentMetadata(1).cigarLength = cigarBuffer_.size() - templ.getFragmentMetadata(1).cigarOffset;
}


