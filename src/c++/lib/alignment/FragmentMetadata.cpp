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
 ** \file FragmentBuilder.cpp
 **
 ** \brief See FragmentMetadata.hh
 ** 
 ** \author Roman Petrovski
 **/

#include "alignment/FragmentMetadata.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{

/**
 * \brief process position adjustment components such as deletions and back dels normally following the contig change
 *
 * @i points to cigar component causing contig change, points before the first normal component on return
 */
int64_t FragmentMetadata::processBackDels(const AlignmentCfg &cfg, const Cigar &cigarBuffer, const unsigned cigarEnd, unsigned &cigarOffset)
{
    int64_t ret = 0;
    for (++cigarOffset; cigarOffset < cigarEnd; ++cigarOffset)
    {
        const std::pair<int, Cigar::OpCode> cigar = Cigar::decode(cigarBuffer[cigarOffset]);
        switch (cigar.second)
        {
            case Cigar::DELETE:
            {
                ret += cigar.first;
                this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (cigar.first - 1) * cfg.normalizedGapExtendScore_);
                break;
            }

            case Cigar::BACK:
            {
                ret -= cigar.first;
                this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (cigar.first - 1) * cfg.normalizedGapExtendScore_);
                break;
            }

            case Cigar::SOFT_CLIP:
            {
                --cigarOffset;
                return ret;
            }

            case Cigar::HARD_CLIP:
            {
                --cigarOffset;
                return ret;
            }

            case Cigar::ALIGN:
            {
                --cigarOffset;
                return ret;
            }

            default:
            {
                ISAAC_ASSERT_MSG(false, "Reached unexpected CIGAR element before reaching ALIGN element." << *this);
            }
        }
    }
    ISAAC_ASSERT_MSG(false, "Reached the end of CIGAR before reaching ALIGN element." << *this);
    return -1UL;
}

//
//unsigned updateAnchors(
//    unsigned sequenceOffset,
//    const unsigned length,
//    reference::Contig::const_iterator currentReference,
//    std::vector<char>::const_iterator currentSequence,
//    Anchor& firstAnchor, Anchor& lastAnchor)
//{
//    static const unsigned MIN_ANCHOR_LENGTH = 32;
//    unsigned ret = 0;
//    unsigned matchesInARow = 0;
//    // scan backwards so that it's easier to check for k-uniqueness
//    currentReference += length;
//    sequenceOffset += length;
//    currentSequence += length;
//    for (unsigned j = 0; length > j; ++j)
//    {
//        --currentReference;
//        --sequenceOffset;
//        --currentSequence;
//        if (isMatch(*currentSequence, *currentReference))
//        {
//            ++matchesInARow;
//            if (matchesInARow >= MIN_ANCHOR_LENGTH)
//            {
//                if (firstAnchor.empty() || firstAnchor.first > sequenceOffset)
//                {
//                    firstAnchor.first = sequenceOffset;
//                    firstAnchor.second = sequenceOffset + MIN_ANCHOR_LENGTH;
//                }
//                if (lastAnchor.empty() || lastAnchor.first < sequenceOffset)
//                {
//                    lastAnchor.first = sequenceOffset;
//                    lastAnchor.second = sequenceOffset + MIN_ANCHOR_LENGTH;
//                }
//            }
//        }
//        else
//        {
//            ret = std::max(ret, matchesInARow);
//            matchesInARow = 0;
//        }
//    }
//    return ret;
//}

double FragmentMetadata::calculateLogProbability(
    unsigned length,
    reference::Contig::const_iterator currentReference,
    std::vector<char>::const_iterator currentSequence,
    std::vector<char>::const_iterator currentQuality)
{
    double ret = 0.0;
    while (length--)
    {
        if (isMatch(*currentSequence, *currentReference))
        {
            ret += Quality::getLogMatch(*currentQuality);
        }
        else
        {
            ret += Quality::getLogMismatch(*currentQuality);
        }
        ++currentReference;
        ++currentSequence;
        ++currentQuality;
    }
    return ret;
}

double FragmentMetadata::calculateInsertionLogProbability(
    unsigned length,
    std::vector<char>::const_iterator currentQuality) const
{
//    // assume the insertion completes the reference by being introduced exactly as in this read
    double ret  = std::accumulate(currentQuality, currentQuality + length, 0.0,
                           bind(std::plus<double>(), _1, boost::bind(&Quality::getLogMatch, _2)));

//    const char qualityMax = *std::max_element(currentQuality, currentQuality + length);
//    // assume one highest-quality base mismatches
//    ret -= Quality::getLogMatch(qualityMax);
//    ret += Quality::getLogMismatch(qualityMax);
    return ret;
}

/**
 * \brief assume a deletion is at least worth a highest base quality mismatch, otherwise one-mismatch alignments
 *        get thrown away for zero-mismatch gapped madness
 */
double calculateDeletionLogProbability(
    std::vector<char>::const_iterator qualityBegin,
    std::vector<char>::const_iterator qualityEnd)
{
//    const char qualityMax = *std::max_element(qualityBegin, qualityEnd);
//    return Quality::getLogMismatch(qualityMax);
    return 0.0;
}

void FragmentMetadata::addMismatchCycles(
    std::vector<char>::const_iterator currentSequence,
    reference::Contig::const_iterator currentReference,
    unsigned sequenceOffset,
    unsigned length, bool reverse,
    const unsigned lastCycle,
    const unsigned firstCycle)
{
    while (length--)
    {
        if (!isMatch(*currentSequence, *currentReference))
        {
            addMismatchCycle(reverse ? lastCycle - sequenceOffset : firstCycle + sequenceOffset);
        }
        ++currentReference;
        ++sequenceOffset;
        ++currentSequence;
    }
}

inline void fixupAnchor(
    const unsigned beginClipOffset,
    const unsigned endClipOffset,
    Anchor& anchor)
{
    anchor.first = std::max<unsigned>(anchor.first, beginClipOffset);
    anchor.second = std::max(anchor.first, anchor.second);
    anchor.second = std::min<unsigned>(anchor.second, endClipOffset);
    anchor.first = std::min(anchor.first, anchor.second);
}

void FragmentMetadata::processInsertion(
    const AlignmentCfg& cfg,
    unsigned &offset,
    const unsigned length,
    std::vector<char>::const_iterator qualityBegin)
{
    this->logProbability += calculateInsertionLogProbability(length, qualityBegin + offset);
    // this discriminates against insertions of mismatching bases making some chimeric alignments look better than proper pairs with insertions in the rescued mate
    //            this->logProbability += calculateLogProbability(
    //                arg, referenceBegin + currentPosition, sequenceBegin + currentBase, qualityBegin + currentBase);
    this->editDistance += length;
    ++this->gapCount;
    this->smithWatermanScore += cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (length - 1) * cfg.normalizedGapExtendScore_);
//    this->smithWatermanScore += cfg.normalizedGapOpenScore_ + (length - 1) * cfg.normalizedGapExtendScore_;
    offset += length;
}

void FragmentMetadata::processDeletion(
    const AlignmentCfg& cfg,
    unsigned offset,
    const unsigned length,
    const std::vector<char>::const_iterator qualityBegin,
    const std::vector<char>::const_iterator qualityEnd,
    int64_t &currentPosition)
{
    this->logProbability += calculateDeletionLogProbability(qualityBegin + offset, qualityEnd);
    this->editDistance += length;
    ++this->gapCount;
    this->smithWatermanScore +=
        cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (length - 1) * cfg.normalizedGapExtendScore_);
    this->splitAlignment |= (length > cfg.splitGapLength_);
    currentPosition += length;
}

void FragmentMetadata::processNegativeDeletion(
    const AlignmentCfg& cfg,
    unsigned offset,
    const unsigned length,
    const std::vector<char>::const_iterator qualityBegin,
    const std::vector<char>::const_iterator qualityEnd,
    int64_t &currentPosition)
{
    this->logProbability += calculateDeletionLogProbability(qualityBegin + offset, qualityEnd);
    this->editDistance += length;
    ++this->gapCount;
    this->smithWatermanScore +=
        cfg.normalizedGapOpenScore_ + std::min(cfg.normalizedMaxGapExtendScore_, (length - 1) * cfg.normalizedGapExtendScore_);
    this->splitAlignment = true;
    currentPosition -= length;
}

void FragmentMetadata::processSoftClip(
    unsigned &offset,
    const unsigned length,
    const std::vector<char>::const_iterator qualityBegin)
{
    // With inversions, soft clipping can occur in the middle of CIGAR
    //            ISAAC_ASSERT_MSG(0 == i || i + 1 == this->cigarLength, "Soft clippings are expected to be "
    //                "found only at the ends of cigar string");
    this->logProbability = std::accumulate(qualityBegin + offset, qualityBegin + offset + length, this->logProbability,
                                           boost::bind(std::plus<double>(), _1, boost::bind(Quality::getLogMatch, _2)));
    // NOTE! Not advancing the reference for soft clips
    offset += length;
}

void FragmentMetadata::processHardClip(
    unsigned &offset,
    const unsigned length)
{
    // Currently hard clips can only be found in inversions. They mask the part of the read
    // that was mapped at the other side of inversion.
    offset += length;
}

void FragmentMetadata::processContigChange(
    const AlignmentCfg& cfg,
    const unsigned newContigId,
    unsigned &currentCigarOffset,
    const reference::ContigList& contigList,
    reference::Contig::const_iterator& referenceBegin,
    unsigned &currentContigId,
    int64_t &currentPosition)
{
    currentPosition += processBackDels(cfg, *this->cigarBuffer, this->cigarOffset + this->cigarLength, currentCigarOffset);
    ISAAC_ASSERT_MSG(0 <= currentPosition, "Unexpected negative position adjustment: " << currentPosition << " " << *this);
    ISAAC_ASSERT_MSG(contigList.at(newContigId).size() >= std::size_t(currentPosition),
                     "Position adjustment outside the contig bounds: " << currentPosition << " " << *this);
    referenceBegin = contigList.at(newContigId).begin();
    currentContigId = newContigId;
    ++this->gapCount;
    this->smithWatermanScore += cfg.normalizedGapOpenScore_;
    this->splitAlignment = true;
}

void FragmentMetadata::processFlip(
    const AlignmentCfg& cfg,
    unsigned & currentBase,
    const unsigned length,
    bool &currentReverse,
    const Read& read,
    std::vector<char>::const_iterator& sequenceBegin,
    std::vector<char>::const_iterator& qualityBegin,
    std::vector<char>::const_iterator& qualityEnd)
{
    currentReverse = !currentReverse;
    sequenceBegin = read.getStrandSequence(currentReverse).begin();
    qualityBegin = read.getStrandQuality(currentReverse).begin();
    qualityEnd = read.getStrandQuality(currentReverse).end();
    currentBase = read.getLength() - currentBase - length;
    this->smithWatermanScore += cfg.normalizedGapOpenScore_;
    this->splitAlignment = true;
    // Notice, this will count flips followed by CONTIG or position adjustment as multiple gaps which is probably not
    // ideal, but so far the gapCount is not being used for anything that requires precise value.
    ++this->gapCount;
}

unsigned iSAAC_PROFILING_NOINLINE FragmentMetadata::processAlign(
    const AlignmentCfg &cfg,
    unsigned &currentBase,
    const unsigned length,
    const std::vector<char>::const_iterator sequenceBegin,
    const std::vector<char>::const_iterator qualityBegin,
    const reference::Contig::const_iterator referenceBegin,
    const bool collectMismatchCycles,
    const flowcell::ReadMetadata &readMetadata,
    const bool currentReverse,
    int64_t &currentPosition)
{
    const unsigned mismatches = alignment::countMismatches(
        sequenceBegin + currentBase, sequenceBegin + currentBase + length, referenceBegin + currentPosition);
    this->logProbability += calculateLogProbability(
        length, referenceBegin + currentPosition, sequenceBegin + currentBase, qualityBegin + currentBase);

    if (collectMismatchCycles)
    {
        addMismatchCycles(
            sequenceBegin + currentBase, referenceBegin + currentPosition,
            currentBase, length, currentReverse, readMetadata.getLastCycle(), readMetadata.getFirstCycle());
    }

//    const unsigned matches =
//        std::inner_product(
//        sequenceBegin + currentBase, sequenceBegin + currentBase + length, referenceBegin + currentPosition,
//        0, std::plus<unsigned>(), &isMatch);
//

    const unsigned matches = length - mismatches;

    this->mismatchCount += mismatches;
    this->smithWatermanScore += cfg.normalizedMismatchScore_ * mismatches;

//    the edit distance includes all mismatches and ambiguous bases (Ns)
//    this->editDistance += std::inner_product(
//        sequenceBegin + currentBase, sequenceBegin + currentBase + length, referenceBegin + currentPosition,
//        0, std::plus<unsigned>(), std::not_equal_to<char>());

    // current isMatch is direct comparison, no need to count twice
    this->editDistance += mismatches;

    currentPosition += length;
    currentBase += length;
    return matches;
}

unsigned iSAAC_PROFILING_NOINLINE FragmentMetadata::updateAlignment(
    const bool collectMismatchCycles,
    const AlignmentCfg &cfg,
    const flowcell::ReadMetadata &readMetadata,
    const reference::ContigList &contigList,
    bool currentReverse,
    unsigned currentContigId,
    const int64_t strandPosition,
    const Cigar &cigarBuffer,
    const unsigned cigarOffset,
    const unsigned cigarLength)
{
    const Read &read = this->getRead();
    std::vector<char>::const_iterator sequenceBegin = read.getStrandSequence(currentReverse).begin();
    std::vector<char>::const_iterator qualityBegin = read.getStrandQuality(currentReverse).begin();
    std::vector<char>::const_iterator qualityEnd = read.getStrandQuality(currentReverse).end();

    ISAAC_ASSERT_MSG(!contigList.at(currentContigId).empty(), "Reference contig was not loaded for " << *this);

    ISAAC_ASSERT_MSG(0 <= strandPosition, "position must be positive for CIGAR update " << *this << " strandPosition:" << strandPosition <<
                     " CIGAR: " << Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarBuffer.end()));
    reference::Contig::const_iterator referenceBegin = contigList.at(currentContigId).begin();

//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(this->getCluster().getId(), " updateFragmentCigar : " <<
//                                           Cigar::toString(cigarBuffer.begin() + cigarOffset, cigarLength ? cigarBuffer.begin() + cigarOffset + cigarLength: cigarBuffer.end()) << " for " << *this);

    this->reverse = currentReverse;
    this->contigId = currentContigId;
    this->cigarBuffer = &cigarBuffer;
    this->cigarOffset = cigarOffset;
    this->cigarLength = cigarLength ? cigarLength : cigarBuffer.size() - this->cigarOffset;
    // adjust cigarOffset and cigarLength
    ISAAC_ASSERT_MSG(cigarBuffer.size() > this->cigarOffset, "Expecting the new cigar is not empty");

    int64_t currentPosition = strandPosition;
    unsigned currentBase = 0;
    unsigned matchCount = 0;
    for (unsigned currentCigarOffset = this->cigarOffset; this->cigarOffset + this->cigarLength > currentCigarOffset; ++currentCigarOffset)
    {
        const Cigar::Component cigar = Cigar::decode(cigarBuffer[currentCigarOffset]);
        if (cigar.second == Cigar::ALIGN)
        {
            matchCount += processAlign(
                cfg, currentBase, cigar.first, sequenceBegin, qualityBegin, referenceBegin, collectMismatchCycles,
                readMetadata, currentReverse, currentPosition);
        }
        else if (cigar.second == Cigar::INSERT)
        {
            processInsertion(cfg, currentBase, cigar.first, qualityBegin);
        }
        else if (cigar.second == Cigar::DELETE)
        {
            processDeletion(cfg, currentBase, cigar.first, qualityBegin, qualityEnd, currentPosition);
        }
        else if (cigar.second == Cigar::BACK)
        {
            processNegativeDeletion(cfg, currentBase, cigar.first, qualityBegin, qualityEnd, currentPosition);
        }
        else if (cigar.second == Cigar::FLIP)
        {
            processFlip(cfg, currentBase, cigar.first, currentReverse, read, sequenceBegin, qualityBegin, qualityEnd);
        }
        else if (cigar.second == Cigar::CONTIG)
        {
            processContigChange(cfg, cigar.first, currentCigarOffset, contigList, referenceBegin, currentContigId, currentPosition);
        }
        else if (cigar.second == Cigar::SOFT_CLIP)
        {
            processSoftClip(currentBase, cigar.first, qualityBegin);
        }
        else if (cigar.second == Cigar::HARD_CLIP)
        {
            processHardClip(currentBase, cigar.first);
        }
        else
        {
            ISAAC_ASSERT_MSG(false, "Unexpected Cigar OpCode:" << cigar.second);
        }
    }

    this->rStrandPos = reference::ReferencePosition(currentContigId, currentPosition, currentReverse);
    this->position = strandPosition;

    const unsigned endClipOffset = getReadLength() - getEndClippedLength();
    // Make sure empty are not crossing into the soft-clipped ends (including quality trimming).
    // This will mess up split read alignment
    fixupAnchor(getBeginClippedLength(), endClipOffset, firstAnchor_);
    fixupAnchor(getBeginClippedLength(), endClipOffset, lastAnchor_);
//    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(this->getCluster().getId(), " updateFragmentCigar : " << *this);

    return matchCount;
}

} // namespace alignment
} // namespace isaac
