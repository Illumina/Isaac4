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

#ifndef iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH
#define iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH

#include <cppunit/extensions/HelperMacros.h>

class TestGapRealigner : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( TestGapRealigner );
    CPPUNIT_TEST( testFull1 );
    CPPUNIT_TEST( testFull2 );
    CPPUNIT_TEST( testFull3 );
    CPPUNIT_TEST( testFull4 );
    CPPUNIT_TEST( testFull5 );
    CPPUNIT_TEST( testFull6 );
    CPPUNIT_TEST( testFull7 );
    CPPUNIT_TEST( testFull8 );
    CPPUNIT_TEST( testFull9 );
    CPPUNIT_TEST( testFull10 );
    CPPUNIT_TEST( testFull11 );
    CPPUNIT_TEST_SUITE_END();
private:

public:
    TestGapRealigner();
    void setUp();
    void tearDown();

    void testFull1();
    void testFull2();
    void testFull3();
    void testFull4();
    void testFull5();
    void testFull6();
    void testFull7();
    void testFull8();
    void testFull9();
    void testFull10();
    void testFull11();
};

#endif // #ifndef iSAAC_ALIGNMENT_TEST_GAP_REALIGNER_HH

