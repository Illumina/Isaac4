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
 ** \file testSplitReadAligner.hh
 **
 ** Tests for split read aligner.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_TEST_SPLIT_READ_ALIGNER_HH
#define iSAAC_ALIGNMENT_TEST_SPLIT_READ_ALIGNER_HH

#include <cppunit/extensions/HelperMacros.h>

#include <string>

#include "alignment/Cigar.hh"
#include "alignment/FragmentMetadata.hh"
#include "flowcell/Layout.hh"
#include "flowcell/ReadMetadata.hh"
#include "alignment/SequencingAdapter.hh"

class TestSplitReadAligner : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestSplitReadAligner );
    CPPUNIT_TEST( testAll );
    CPPUNIT_TEST_SUITE_END();
private:
    isaac::alignment::Cigar cigarBuffer_;

public:
    TestSplitReadAligner();
    void setUp();
    void tearDown();
    void testAll()
    {
//        ISAAC_SCOPE_BLOCK_CERR
        {
            testEverything1();
            testEverything2();
            testEverything3();
            testEverything4();
            testEverything5();
            testEverything6();
            testEverything7();
            testEverything8();
        }
    }
    void testEverything1();
    void testEverything2();
    void testEverything3();
    void testEverything4();
    void testEverything5();
    void testEverything6();
    void testEverything7();
    void testEverything8();

private:
    void align(
        const std::string &read,
        const std::string &reference,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList,
        const unsigned clusterId = 0);

    void align(
        const std::string &readAlignment1,
        const std::string &readAlignment2,
        const std::string &reference,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList,
        const unsigned clusterId = 0);

    void align(
        const std::string &readAlignment1,
        const std::string &readAlignment2,
        const std::string &reference,
        const std::string &reference2,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList,
        const unsigned clusterId = 0);

    void align(
        const isaac::alignment::Cluster &cluster,
        const isaac::flowcell::ReadMetadataList &readMetadatList,
        const std::vector<char> reference1WithoutSpaces,
        const std::vector<char> reference2WithoutSpaces,
        isaac::alignment::FragmentMetadataList &fragmentMetadataList);
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_SPLIT_READ_ALIGNER_HH

