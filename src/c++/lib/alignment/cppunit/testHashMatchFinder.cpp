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
 ** \file TestHashMatchFinder.cpp
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
#include <boost/bind.hpp>
#include <boost/lambda/construct.hpp>

using namespace std;

#include "RegistryName.hh"
#include "testHashMatchFinder.hh"

#include "alignment/HashMatchFinder.hh"
#include "alignment/matchFinder/TileClusterInfo.hh"
#include "BuilderInit.hh"
#include "common/Threads.hpp"
#include "reference/SortedReferenceMetadata.hh"
#include "reference/ReferenceHasher.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestHashMatchFinder, registryName("HashMatchFinder"));

TestHashMatchFinder::TestHashMatchFinder()
{

}

void TestHashMatchFinder::setUp()
{
}

void TestHashMatchFinder::tearDown()
{
}

TestMatchStorage TestHashMatchFinder::findMatches(
    const std::string& reference, const std::string& sequence,
    const isaac::flowcell::ReadMetadataList &readMetadataList)
{
    TestContigList contigList(reference);

    isaac::common::ThreadVector threads(1);
    isaac::reference::ReferenceHasher<isaac::reference::ReferenceHash<isaac::oligo::VeryShortKmerType> > referenceHasher(
        contigList, threads, threads.size());

    const isaac::reference::ReferenceHash<isaac::oligo::VeryShortKmerType> referenceHash = referenceHasher.generate(0x10000);

    isaac::flowcell::FlowcellLayoutList flowcells(1, isaac::flowcell::Layout("", isaac::flowcell::Layout::Fastq, isaac::flowcell::FastqFlowcellData(false, '!', false), 8, 0, std::vector<unsigned>(),
                                         readMetadataList, "blah"));

    isaac::flowcell::TileMetadata tileMetadata(
        flowcells.front().getFlowcellId(), flowcells.front().getIndex(), 0, 1, 1, 0);
    const isaac::flowcell::TileMetadataList tileMetadataList(std::vector<isaac::flowcell::TileMetadata>(1, tileMetadata));

    const unsigned clusterLength = isaac::flowcell::getTotalReadLength(readMetadataList);
    isaac::alignment::BclClusters tileClusters(clusterLength);

    tileClusters.reset(clusterLength, 1);

//    const std::vector<char>& bcls = getBcl(readMetadataList, contigList, 0, 0, readMetadataList.front().getLength());
    const std::vector<char>& bcls = getBcl(sequence.substr(0, clusterLength));

    std::copy(bcls.begin(), bcls.end(), tileClusters.cluster(0));

    isaac::alignment::Cluster cluster(readMetadataList.front().getLength());
    cluster.init(flowcells.front().getReadMetadataList(), tileClusters.cluster(0),
                 0, 0, isaac::alignment::ClusterXy(0,0), true, 0, 0);


    TestMatchStorage matchLists(100);
//    isaac::alignment::ClusterHashMatchFinder<isaac::oligo::VeryShortKmerType> matchFinder(
//        referenceHash, flowcells, isaac::flowcell::BarcodeMetadataList(), 0, repeatThreshold, std::vector<std::size_t>(),
//        sortedReferenceMetadataList, seedMetadataList, seedMetadataList.size());
    isaac::alignment::ClusterHashMatchFinder< isaac::reference::ReferenceHash<isaac::oligo::VeryShortKmerType>, 4> matchFinder(
        referenceHash, 1000, 0, 1000);

    isaac::alignment::ReferenceOffsetLists fwMergeBuffers(11, isaac::alignment::ReferenceOffsetList(1000));
    isaac::alignment::ReferenceOffsetLists rvMergeBuffers(11, isaac::alignment::ReferenceOffsetList(1000));

    isaac::alignment::matchFinder::TileClusterInfo tileClusterInfo(tileMetadataList);
    tileClusterInfo.setBarcodeIndex(0, 0, 0);
    matchFinder.findReadMatches(
        contigList, cluster, flowcells.front().getReadMetadataList().front(), 1000, matchLists, fwMergeBuffers, rvMergeBuffers);

    return matchLists;
}

void TestHashMatchFinder::testEverything()
{
//    ISAAC_SCOPE_BLOCK_CERR
    {

    {
        std::string reference("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        isaac::flowcell::ReadMetadataList readMetadataList(1, isaac::flowcell::ReadMetadata(1, sequence.length() + 1, 0, 0));
        TestMatchStorage matchLists = findMatches(reference, sequence, readMetadataList);

        CPPUNIT_ASSERT_EQUAL(1UL, matchLists.at(3).size());
    }

    {
        std::string reference("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAA"
                              "ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        isaac::flowcell::ReadMetadataList readMetadataList(1, isaac::flowcell::ReadMetadata(1, sequence.length() + 1, 0, 0));
        TestMatchStorage matchLists = findMatches(reference, sequence, readMetadataList);

        CPPUNIT_ASSERT_EQUAL(1UL, matchLists.at(3).size());
    }

    {
        std::string reference("TTTTTTGTTTTTTCTTCTTGTTATCTTTATTTTTTTAAT"
                              "ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
        isaac::flowcell::ReadMetadataList readMetadataList(1, isaac::flowcell::ReadMetadata(1, sequence.length() + 1, 0, 0));
        TestMatchStorage matchLists = findMatches(reference, sequence, readMetadataList);

        CPPUNIT_ASSERT_EQUAL(1UL, matchLists.at(3).size());
    }

    {
        std::string reference("GTGGGGGAAGCTGAGTCTCACTTTGTCGCCCAGGCTGGAGTGCAGCGGCGCCATTTCAGCTCACTGTAACCTCCACCTCTGTGATTCAAGCAATTCTCAT");
        std::string sequence ("GTGGGGGAAGCTGAGTCTCACTTTGTCGCCCAGGCTGGAGTGCAGCGGCGCCATTTCAGCTCACTGTAACCTCCACCTCTGTGATTCAAGCAATTCTCAT");
        isaac::flowcell::ReadMetadataList readMetadataList(1, isaac::flowcell::ReadMetadata(1, sequence.length() + 1, 0, 0));
        TestMatchStorage matchLists = findMatches(reference, sequence, readMetadataList);

        CPPUNIT_ASSERT_EQUAL(0UL, matchLists.at(1).size());
        CPPUNIT_ASSERT_EQUAL(1UL, matchLists.at(3).size());
        CPPUNIT_ASSERT_EQUAL(1000U, matchLists.at(3).at(0).contigListOffset_);
        CPPUNIT_ASSERT_EQUAL(false, matchLists.at(3).at(0).reverse_);
    }

//    current iterative match finder will not stop until it finds 3 seeds support. 16 bases is not enough
//    //C C G T C C C C
//    {
//        std::cerr << "tada" << std::endl;
//        std::string reference("ACTCAGCTTCCCCCAC");
//        std::string sequence ("GTGGGGGAAGCTGAGT");
//        isaac::flowcell::ReadMetadataList readMetadataList(1, isaac::flowcell::ReadMetadata(1, 16, 0, 0));
//        TestMatchStorage matchLists = findMatches(reference, sequence, readMetadataList);
//
//        CPPUNIT_ASSERT_EQUAL(4UL, matchLists.size());
//        CPPUNIT_ASSERT_EQUAL(0UL, matchLists.at(1).size());
//        CPPUNIT_ASSERT_EQUAL(1UL, matchLists.at(3).size());
//        CPPUNIT_ASSERT_EQUAL(1000U, matchLists.at(3).at(0).contigListOffset_);
//        CPPUNIT_ASSERT_EQUAL(true, matchLists.at(3).at(0).reverse_);
//    }

    //C C G T C C C C
    {
        std::string reference("ATGAGAATTGCTTGAATCACAGAGGTGGAGGTTACAGTGAGCTGAAATGGCGCCGCTGCACTCCAGCCTGGGCGACAAAGTGAGACTCAGCTTCCCCCAC");
        std::string sequence ("GTGGGGGAAGCTGAGTCTCACTTTGTCGCCCAGGCTGGAGTGCAGCGGCGCCATTTCAGCTCACTGTAACCTCCACCTCTGTGATTCAAGCAATTCTCAT");
        isaac::flowcell::ReadMetadataList readMetadataList(1, isaac::flowcell::ReadMetadata(1, 100, 0, 0));
        TestMatchStorage matchLists = findMatches(reference, sequence, readMetadataList);

        CPPUNIT_ASSERT_EQUAL(0UL, matchLists.at(1).size());
        CPPUNIT_ASSERT_EQUAL(1UL, matchLists.at(3).size());
        CPPUNIT_ASSERT_EQUAL(1000U, matchLists.at(3).at(0).contigListOffset_);
        CPPUNIT_ASSERT_EQUAL(true, matchLists.at(3).at(0).reverse_);
    }
// SeedMetadataList is not used anymore.
//    {
//        std::string reference("TTTTTTGTTTTTTCTTCTTGTTATCTTTATTTTTTTAAT"
//                              "ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 10,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);
//
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
//    }

// SeedMetadataList is not used anymore.
//    {
//        std::string reference(//"TTTTTTGTTTTTTCTTCTTGTTATCTTTATTTTTTTAAT"
//                              "ATTAAAAAAATAAAGA"
//                              "ATTAAAAAAATAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 2,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);
//
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), !boost::bind(&isaac::alignment::Match::isTooManyMatch, _1))));
//    }

// SeedMetadataList is not used anymore.
//    {
//        std::string reference("TTTTTTGTATTTTCTTCTTGTTA"
//            // reverse
//                              "TCTTTACATTTTTAAT"
//            // forward:
//                              "ATTAAAAATGTAAAGT"
//                              "TAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        std::string sequence ("ATTAAAAATGTAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 2,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);
//
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
//
////        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 39, false, true), matches.at(1).location_);
////        CPPUNIT_ASSERT_EQUAL(1U, matches.at(1).matchCount_);
////
////        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 89, false, false), matches.at(0).location_);
////        CPPUNIT_ASSERT_EQUAL(1U, matches.at(0).matchCount_);
//    }

// SeedMetadataList is not used anymore.
//    {
//        std::string reference("TTTTTTGTATTTTCTTCTTGTTA"
//            // reverse
//                              "TCTTTACATTTTTAAT"
//            // forward:
//                              "ATTAAAAATGTAAAGT"
//                              "TAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        std::string sequence ("ATTAAAAATGTAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 2,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0))(isaac::alignment::SeedMetadata( 4, 8, 0, 1)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);
//
//        for (const auto &m: matches)
//        {
//            ISAAC_THREAD_CERR << "blah: " << m << std::endl;
//        }
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
//
//        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 89, false, false), matches.at(0).location_);
//        // seeds overlap. So, one match cound only
//        CPPUNIT_ASSERT_EQUAL(1U, matches.at(0).seedCount_);
//
//        CPPUNIT_ASSERT_EQUAL(isaac::reference::ReferencePosition(0, 39, false, true), matches.at(1).location_);
//        // seeds overlap. So, one match cound only
//        CPPUNIT_ASSERT_EQUAL(1U, matches.at(1).seedCount_);
//    }

// SeedMetadataList is not used anymore.
//    {
//        std::string reference("GGGGTTTTCCCAACCTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        std::string sequence ("ATTAAAAAAATAAAGATAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 2,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( sequence.length() - 8, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(),
//            boost::assign::list_of(isaac::flowcell::ReadMetadata( 1, sequence.length() + 1, 0, 0)).convert_to_container<isaac::flowcell::ReadMetadataList>());
//
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
//    }

// SeedMetadataList is not used anymore.
//    {
//        std::string reference("ATTAAAAAATTAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTATTAAAAA");
//        std::string sequence ("ATTAAAAAATTAAAGTTAACAAGAAGAAAAAACAAAAAACAGAAAATAATTAAACAGGGACAAACCAAAGACAAAATACGATTTGGAAGAAGGCCACAAAAAAACCCCTTTAGGGGGGGTTTTCCCAACC");
//        TestMatchStorage matches = findMatches(
//            reference, sequence, 2,
//            boost::assign::list_of(isaac::alignment::SeedMetadata( 0, 8, 0, 0)).convert_to_container<isaac::alignment::SeedMetadataList>(), readMetadataList_);
//
//        CPPUNIT_ASSERT_EQUAL(2UL, std::size_t(std::count_if(matches.begin(), matches.end(), boost::bind(&isProperMatch, _1))));
//    }
    }
}
