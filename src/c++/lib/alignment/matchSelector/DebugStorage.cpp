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
 ** \file DebugFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <stdio.h>
#include "alignment/matchSelector/DebugStorage.hh"

#include "demultiplexing/BarcodePathMap.hh"
#include "alignment/matchSelector/debugStorage/Utils.hh"
#include "common/Debug.hh"
#include "common/FastIo.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{

const int READ_LENGTH_MAX = 1000;

DebugStorage::DebugStorage(
    const reference::ContigList &contigList,
    const AlignmentCfg &alignmentCfg,
    const flowcell::FlowcellLayoutList &flowcellLayoutList,
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const unsigned threads,
    alignment::matchSelector::FragmentStorage &actualStorage):
    contigList_(contigList),
    flowcell_(flowcellLayoutList.at(0)),
    alignmentCfg_(alignmentCfg),
    outputDirectory_(outputDirectory),
    barcodeQQR1Paths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "R1QQ.tsv")),
    barcodeQQR2Paths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "R2QQ.tsv")),
    barcodeQQR1OriPaths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "R1QQOri.tsv")),
    barcodeQQR2OriPaths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "R2QQOri.tsv")),
    barcodeMapqPaths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "Mapq-00.tsv")),
    barcodeSMPaths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "SM-00.tsv")),
    barcodeASPaths_(demultiplexing::mapBarcodesToFiles(outputDirectory_, barcodeMetadataList, "AS-00.tsv")),
    actualStorage_(actualStorage),
    originalCigars_(READ_LENGTH_MAX)
{
    createDirectories(barcodeQQR1Paths_, barcodeMetadataList);

//    while (pairMapqStatistics_.size() < classesCount)
//    {
//        char szDebugClass[100];
//        sprintf(szDebugClass, "%02x", unsigned(pairMapqStatistics_.size()));
//        pairMapqStatistics_.push_back(new debugStorage::MapqStatistics(
//            pairMapqStatistics_.size(),
//            outputDirectory_ / ("MapQ-" + std::string(szDebugClass) + ".tsv"),
//            outputDirectory_ / ("SM-" + std::string(szDebugClass) + ".tsv"),
//            outputDirectory_ / ("AS-" + std::string(szDebugClass) + ".tsv")));
//
//        r1QqOriginalStatistics_(contigList_, outputDirectory / "R1QQOri.tsv"),
//        r2QqOriginalStatistics_(contigList_, outputDirectory / "R2QQOri.tsv"),
//        r1QqStatistics_(contigList_, outputDirectory / "R1QQ.tsv"),
//        r2QqStatistics_(contigList_, outputDirectory / "R2QQ.tsv"),
//    }

    for (unsigned sampleId = 0; sampleId < barcodeQQR1Paths_.getTotalSamples(); ++sampleId)
    {
        pairMapqStatistics_.push_back(new debugStorage::MapqStatistics(
                pairMapqStatistics_.size(),
                barcodeMapqPaths_.getSampleFilePath(sampleId),
                barcodeSMPaths_.getSampleFilePath(sampleId),
                barcodeASPaths_.getSampleFilePath(sampleId)));
        r1QqStatistics_.push_back(new debugStorage::QqStatistics(contigList_, barcodeQQR1Paths_.getSampleFilePath(sampleId)));
        r2QqStatistics_.push_back(new debugStorage::QqStatistics(contigList_, barcodeQQR2Paths_.getSampleFilePath(sampleId)));
        r1QqOriginalStatistics_.push_back(new debugStorage::QqStatistics(contigList_, barcodeQQR1OriPaths_.getSampleFilePath(sampleId)));
        r2QqOriginalStatistics_.push_back(new debugStorage::QqStatistics(contigList_, barcodeQQR2OriPaths_.getSampleFilePath(sampleId)));
    }

    originalCigars_.at(0).resize(threads);
    originalCigars_.at(1).resize(threads);
    // make all sorts of cigars so that we don't need to allocate memory during runtime and keep memory control happy
    for (unsigned t = 0; t < threads; ++t)
    {
        originalCigars_.at(0).at(t).reserve(100);
        originalCigars_.at(1).at(t).reserve(100);
    }
}

DebugStorage::~DebugStorage()
{
}

bool DebugStorage::restoreOriginal(
    const unsigned threadNumber,
    const std::size_t readNumber,
    FragmentMetadata &fragment)
{
    const reference::ReferencePosition oriPos = debugStorage::getAlignmentPositionFromName(readNumber, fragment);
    if (oriPos.isTooManyMatch())
    {
        return false;
    }

    std::pair<const char *, const char *> nameCigar = debugStorage::getCigarFromName(readNumber, fragment);
    originalCigars_.at(readNumber - 1).at(threadNumber).fromString(nameCigar.first, nameCigar.second);


    fragment.resetAlignment();
    if (!oriPos.isNoMatch())
    {
        fragment.updateAlignment(
            false,
            alignmentCfg_,
            flowcell_.getReadMetadataList().at(fragment.getReadIndex()),
            contigList_,
            oriPos.reverse(),
            oriPos.getContigId(),
            oriPos.getPosition(),
            originalCigars_.at(readNumber - 1).at(threadNumber),
            0);
    }
    return true;
}


void DebugStorage::updateMapqStats(
    const BamTemplate& bamTemplate,
    const unsigned barcodeIdx)
{
    const isaac::alignment::FragmentMetadata &r0 = bamTemplate.getFragmentMetadata(0);
    const unsigned int firstReadNumber = flowcell_.getReadMetadataList().at(r0.getReadIndex()).getNumber();
    if (2 == bamTemplate.getFragmentCount())
    {
        const isaac::alignment::FragmentMetadata &r1 = bamTemplate.getFragmentMetadata(1);
        unsigned int secondReadNumber = flowcell_.getReadMetadataList().at(r1.getReadIndex()).getNumber();

        pairMapqStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(bamTemplate, firstReadNumber, r0, &r1);
        pairMapqStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(bamTemplate, secondReadNumber, r1, &r0);
    }
    else
    {
        pairMapqStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(bamTemplate, firstReadNumber, r0, 0);
    }
}


void DebugStorage::store(
    const BamTemplate &bamTemplate,
    const unsigned barcodeIdx,
    const unsigned threadNumber)
{
    updateMapqStats(bamTemplate, barcodeIdx);
    BamTemplate originalTemplate = bamTemplate;
    bool originalRestored = true;
    const isaac::alignment::FragmentMetadata &r1Aligned = bamTemplate.getFragmentMetadata(0);
    const unsigned int r1Number = flowcell_.getReadMetadataList().at(r1Aligned.getReadIndex()).getNumber();
    bool changed = false;//!r1.isAligned();

//    ISAAC_THREAD_CERR<< "Checking: " << r1Aligned << " " << originalTemplate.getAlignmentScore() << "AS" << std::endl;

    if (r1Aligned.hasMapQ() &&
        !alignment::containsHomopolymer(r1Aligned.getStrandSequence().begin(), r1Aligned.getStrandSequence().end()) &&
        !debugStorage::alignsCorrectly(r1Number, r1Aligned))
    {
        changed = true;
        if (r1Aligned.mapQ >= 40)
        {
            ISAAC_THREAD_CERR<< "Misaligned: " << r1Aligned << " " << originalTemplate.getAlignmentScore() << "AS" << std::endl;
        }
    }

    isaac::alignment::FragmentMetadata &r1 = originalTemplate.getFragmentMetadata(0);
    r1QqStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(0, bamTemplate, true);
    if (restoreOriginal(threadNumber, flowcell_.getReadMetadataList().at(r1.getReadIndex()).getNumber(), r1))
    {
        r1QqOriginalStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(0, originalTemplate, false);
    }
    else
    {
        originalRestored = false;
    }
    if (2 == bamTemplate.getFragmentCount())
    {
        isaac::alignment::FragmentMetadata &r2 = originalTemplate.getFragmentMetadata(1);
        r2QqStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(1, bamTemplate, true);
        // use aligned, not original r1
        const isaac::alignment::FragmentMetadata &r2Aligned = bamTemplate.getFragmentMetadata(1);
        unsigned int r2Number = flowcell_.getReadMetadataList().at(r2Aligned.getReadIndex()).getNumber();
//        ISAAC_THREAD_CERR<< "Checking: " << r2Aligned << " " << originalTemplate.getAlignmentScore() << "AS" << std::endl;
        if (r2Aligned.hasMapQ() &&
            !alignment::containsHomopolymer(r2Aligned.getStrandSequence().begin(), r2Aligned.getStrandSequence().end()) &&
            !debugStorage::alignsCorrectly(r2Number, r2Aligned))
        {
            changed = true;
            if (r2Aligned.mapQ >= 40)
            {
                ISAAC_THREAD_CERR<< "Misaligned: " << r2Aligned << " " << originalTemplate.getAlignmentScore() << "AS" << std::endl;
            }
        }

//        if (!r1Aligned.isAligned() && !r2Aligned.isAligned())
//        {
//            ISAAC_THREAD_CERR<< "Unaligned: " << r1 << " " << r2 << std::endl;
//        }
//        store = r1Aligned.isAligned() != r2.isAligned();// && debugStorage::alignsCorrectly(r1Number, r1Aligned) && !debugStorage::alignsCorrectly(r2Number, r2);
//        store = r1Aligned.isAligned() && r2Aligned.isAligned() && (
//            (!debugStorage::alignsCorrectly(r1Number, r1Aligned) && r1Aligned.mapQ >= 30) ||
//            (!debugStorage::alignsCorrectly(r2Number, r2Aligned) && r2Aligned.mapQ >= 30));
//        store =
//            (r1Aligned.isAligned() && !debugStorage::alignsCorrectly(r1Number, r1Aligned) && 46 <= r1Aligned.mapQ && r1Aligned.mapQ <= 58) ||
//            (r2.isAligned() && !debugStorage::alignsCorrectly(r2Number, r2) && 46 <= r2.mapQ && r2.mapQ <= 58);
//        store = store || !r2.isAligned();

//        store = r1Aligned.isAligned() && debugStorage::alignsCorrectly(r1Number, r1Aligned) && r2Aligned.isAligned() && debugStorage::alignsCorrectly(r2Number, r2Aligned) &&
//            45 == r1Aligned.mapQ;

//        store |= r2.getCluster().getId() == 1761956 || r2Aligned.gapCount || (r2Aligned.isAligned() && debugStorage::alignsCorrectly(r2Number, r2Aligned) && r2Aligned.mismatchCount > 20 && r2Aligned.secondBestMismatchDelta);

        if (restoreOriginal(threadNumber, r2Number, r2))
        {
            r2QqOriginalStatistics_[barcodeQQR1Paths_.getSampleIndex(barcodeIdx)].updateStat(1, originalTemplate, false);
        }
        else
        {
            originalRestored = false;
        }
    }

//    store = r1Aligned.isAligned() && debugStorage::alignsCorrectly(r1Number, r1Aligned) && 45 <= r1Aligned.mapQ && r1Aligned.mapQ <= 59;
//    store |= r1.getCluster().getId() == 1761956 || r1Aligned.gapCount || (r1Aligned.isAligned() && debugStorage::alignsCorrectly(r1Number, r1Aligned) && r1Aligned.mismatchCount > 20 && r1Aligned.secondBestMismatchDelta);
    actualStorage_.store(bamTemplate, barcodeIdx, threadNumber);
    if (changed)
    {
//        actualStorage_.store(originalTemplate, barcodeIdx);
        actualStorage_.store(bamTemplate, barcodeIdx + 2, threadNumber);
        if (originalRestored)
        {
            actualStorage_.store(originalTemplate, barcodeIdx + 3, threadNumber);
        }
    }
    else
    {
        // store unchanged
        actualStorage_.store(bamTemplate, barcodeIdx + 1, threadNumber);
    }
}

} //namespace matchSelector
} // namespace alignment
} // namespace isaac
