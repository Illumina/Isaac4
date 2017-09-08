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
 ** \file FragmentMetadata.hh
 **
 ** \brief Component to encapsulate all metadata associated to a fragment.
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH
#define iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH

#include <vector>
#include <algorithm>
#include <iostream>
#include <bitset>

#include "../common/StaticVector.hh"
#include "alignment/AlignmentCfg.hh"
#include "alignment/Mismatch.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/Quality.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace alignment
{


struct Anchor : public std::pair<unsigned short, unsigned short>
{
    Anchor(unsigned short f, unsigned short s, bool kUnique) : std::pair<unsigned short, unsigned short>(f,s), kUnique_(kUnique)//, kRepeat_(false)
        {}
    Anchor(bool kUnique = false) : std::pair<unsigned short, unsigned short>(0,0), kUnique_(kUnique)//, kRepeat_(false)
        {}
    unsigned short length() const {return second - first;}
    bool empty() const {return second == first;}
    bool good() const {return !empty()/* && (kUnique_ || kRepeat_)*/;}
    bool betterThan(const Anchor &that) const {return (!empty() && that.empty()) || (good() && !that.good());}

    bool kUnique_ : 1;
//    bool kRepeat_ : 1;
//
    friend std::ostream &operator <<(std::ostream &os, const Anchor &a)
    {
        return os << "Anchor(" << a.first << "," << a.second
             << "," << (a.kUnique_ ? "ku," : "")
            // << (a.kRepeat_ ? "kr" : "")
            << ")";
    }
};

/**
 ** \brief Alignment information for a fragment (as defined by the SAM Format
 ** Specification (v1.4-r962)
 **
 ** This component is the building block for the FragmentBuilder. It is
 ** designed for efficiency, and does not involve any memory allocation (which
 ** is the reason why it does not and should not store the CIGAR).
 **/
struct FragmentMetadata
{
    static const unsigned MAX_CYCLES = 1024;
    static const unsigned POSSIBLY_SEMIALIGNED_MISMATCH_COUNT = 8;
    /**
     ** \brief Cluster associated to the fragment
     **
     ** Note, it is the responsibility of the calling code to ensure the
     ** validity of the cluster.
     **/
    const Cluster *cluster;

    /// Id of the contig where the fragment is located
    unsigned contigId;
    /**
     ** \brief 0-based leftmost position of the fragment on the forward strand
     ** of the contig.
     **
     ** Even though the position can become negative while building the fragment
     ** (before calculating the cigar), the final position will be guaranteed to
     ** be positive (or 0). The final position is the position of the first
     ** aligned base ('M', '=' or 'X'). If the read extends outside the contig,
     ** this will be reflected by appropriate insertions or clipping operations
     ** in the cigar.
     **/
    int64_t position;

    /// number of bases clipped from the lowest read cycle irrespective of alignment
    unsigned short lowClipped;
    /// number of bases clipped from the highest read cycle irrespective of alignment
    unsigned short highClipped;

    // number of bases clipped due to adapter
    unsigned short adapterClipped_;

    /**
     ** \brief observed length of the fragment on the contig
     **
     ** If there are no indels and no clipping, this would be the length of the
     ** read. With indels this is the read length minus the insertions plus the
     ** deletions (resp. to and from the reference). Clipped areas of the read
     ** are also subtracted from the length of the read.
     **/
//    unsigned observedLength;
    reference::ReferencePosition rStrandPos;

    /// 0-based index of the read in the list of ReadMetadata
    unsigned readIndex;
    /// Orientation of the read. False is forward, true is reverse
    bool reverse:1;
    bool splitAlignment:1;
    bool decoyAlignment:1;
    // number of equivalent alignments when discovered by seeds if this one is best. 1 means this is best unique alignmet. Otherwise undetermined
    unsigned repeatCount;
    /// Cigar offset in the associated cigar buffer (see FragmentBuilder)
    unsigned cigarOffset;
    /// Number of operations in the cigar
    unsigned cigarLength;
    /// Buffer containing the cigar data
    const Cigar *cigarBuffer;

    /// Number of mismatches in the alignment (can't be more than read length)
    unsigned mismatchCount;

    unsigned uncheckedSeeds;

    /// Number of breakpoints (indels and other splits) in the fragment
    unsigned gapCount;

    /// Edit distance from the alignment (including indels and ambiguous bases)
    unsigned editDistance;

    bool isSplit() const {return splitAlignment;}

    // set fo cycles containing mismatches (outside indels).
    typedef std::bitset<MAX_CYCLES> MismatchCycles;
    std::bitset<MAX_CYCLES> mismatchCycles;

    class ConstMismatchCycleIterator :
        public boost::iterator_facade<
        ConstMismatchCycleIterator
        , unsigned
        , boost::forward_traversal_tag
        , unsigned const &
        >
    {
        unsigned cycle_;
        boost::reference_wrapper<const MismatchCycles> mismatchCycles_;
     public:
        explicit ConstMismatchCycleIterator(const MismatchCycles &mismatchCycles, const std::size_t cycle = 1)
          : cycle_(cycle), mismatchCycles_(mismatchCycles)
        {
            ISAAC_ASSERT_MSG(cycle_, "Cycle must be a 1-based integer.");
            while (cycle_ != (mismatchCycles_.get().size() + 1) && !mismatchCycles_.get().test(cycle_ - 1)) ++cycle_;
        }

     private:
        friend class boost::iterator_core_access;

        void increment() { ++cycle_; while (cycle_ != (mismatchCycles_.get().size() + 1) && !mismatchCycles_.get().test(cycle_ - 1)) ++cycle_; }

        bool equal(ConstMismatchCycleIterator const& other) const
        {
            ISAAC_ASSERT_MSG(mismatchCycles_.get_pointer() == other.mismatchCycles_.get_pointer(), "Illegal compare for iterators of different containers.");
            return this->cycle_ == other.cycle_;
        }

        const unsigned &dereference() const {
            return cycle_; }
    };

    /**
     ** \brief natural logarithm of the probability that the fragment is correct
     **
     ** This is intrinsic to the fragment and depends only of the quality of all
     ** matching base and and the quality of all mismatching bases (indel are
     ** counted as matching bases). See AlignmentQuality for the detail of the
     ** probabilities used in each case.
     **/
    double logProbability;

    /// lowest subsequence in the direction of the reference that anchors the alignment
    Anchor firstAnchor_;
    /// highest subsequence in the direction of the reference that anchors the alignment
    Anchor lastAnchor_;

    const Anchor &headAnchor() const {return reverse ? lastAnchor_ : firstAnchor_;}
    const Anchor &tailAnchor() const {return reverse ? firstAnchor_ : lastAnchor_;}
    Anchor &headAnchor() {return reverse ? lastAnchor_ : firstAnchor_;}
    Anchor &tailAnchor() {return reverse ? firstAnchor_ : lastAnchor_;}
    /// recompute anchor
    static const unsigned ALLOWED_MISMATCHES_IN_ANCHOR = 3;
    static const unsigned ANCHOR_LENGTH = 16;
    template<bool firstOnly = true>
    bool recomputeTailAnchor(const reference::ContigList& contigList)
        {return recomputeAnchor<false, firstOnly, ALLOWED_MISMATCHES_IN_ANCHOR>(contigList, ANCHOR_LENGTH); }
    template<bool firstOnly = true, unsigned anchorMismatches = ALLOWED_MISMATCHES_IN_ANCHOR>
    bool recomputeHeadAnchor(const reference::ContigList& contigList, const unsigned anchorLength = ANCHOR_LENGTH)
        {return recomputeAnchor<true, firstOnly, anchorMismatches>(contigList, anchorLength);}

    template <bool head,  bool firstOnly, unsigned anchorMismatches>
    const Anchor computeAnchor(const reference::ContigList& contigList, const unsigned anchorLength) const
    {
        ISAAC_ASSERT_MSG(!splitAlignment, "Can't compute anchor on a split alignment:" << *this);
        const Sequence & strandSequence = getStrandSequence();
        const reference::Contig &contig = contigList.at(getContigId());
        if (head != reverse)
        {
            return makeForwardAnchor<firstOnly>(
                contig, anchorLength, anchorMismatches, getPosition(),
                strandSequence.begin(), getBeginClippedLength(), strandSequence.end() - getEndClippedLength());
        }
        else
        {
            return makeReverseAnchor<firstOnly>(
                contig, anchorLength, anchorMismatches, rStrandPos.getPosition(),
                strandSequence, getEndClippedLength(), strandSequence.rend() - getBeginClippedLength());
        }
    }


    /**
     ** \brief Alignment score in the global context of the reference
     **
     ** This depends on the all the pLogCorrect values for all the possible
     ** alignments for this read across the whole reference. It also takes into
     ** account the rest-of-genome correction. Value of -1U indicates unknown alignment score.
     **/
    unsigned alignmentScore;
    unsigned char mapQ;

    // Weighted sum of mismatch and gap penalties similar to what's used for Smith-Waterman alignment
    int smithWatermanScore;

    bool isBetterUngapped(const FragmentMetadata &right) const
    {
        return ISAAC_LP_LESS(right.logProbability, logProbability) ||
            (ISAAC_LP_EQUALS(right.logProbability, logProbability) && smithWatermanScore < right.smithWatermanScore);
    }

    /**
     * \brief prefers small gaps to big or structural variants
     */
    bool isBetterGapped(const FragmentMetadata &right) const
    {
        if (smithWatermanScore < right.smithWatermanScore)
        {
            return true;
        }

        if (smithWatermanScore == right.smithWatermanScore)
        {
            if (splitAlignment < right.splitAlignment)
            {
                return true;
            }
            else if (splitAlignment > right.splitAlignment)
            {
                return false;
            }

            // if both don't have splits pick the shortest.
            if (!splitAlignment && !right.splitAlignment)
            {
                if (getObservedLength() < right.getObservedLength())
                {
                    return true;
                }
                else if (right.getObservedLength() < getObservedLength())
                {
                    return false;
                }
            }
            //split alignment or observed lengths match
            return ISAAC_LP_LESS(right.logProbability, logProbability);
        }

        return false;
    }

    static bool bestUngappedLess(const FragmentMetadata &left, const FragmentMetadata &right)
    {
        return left.isBetterUngapped(right);
    }

    static bool bestGappedLess(const FragmentMetadata &left, const FragmentMetadata &right)
    {
        return left.isBetterGapped(right);
    }

    static bool alignmentsEquivalent(const FragmentMetadata &left, const FragmentMetadata &right)
    {
        return left.smithWatermanScore == right.smithWatermanScore && ISAAC_LP_EQUALS(right.logProbability, left.logProbability);
    }

    /**
     ** \brief Comparison of FragmentMetadata by reference position
     **/
    bool operator<(const FragmentMetadata &f) const
    {
        ISAAC_ASSERT_MSG(cluster == f.cluster && readIndex == f.readIndex,
                         "Comparison makes sense only for metadata representing the same fragment " <<
                         *this << " vs " << f);
        return
            (this->contigId < f.contigId ||
                (this->contigId == f.contigId && (this->position < f.position ||
                    (this->position == f.position && (this->reverse < f.reverse ||
                        (reverse == f.reverse && (rStrandPos < f.rStrandPos ||
                            (rStrandPos == f.rStrandPos &&
                                std::lexicographical_compare(cigarBegin(), cigarEnd(), f.cigarBegin(), f.cigarEnd())
                            ))))))));
    }

    bool operator == (const FragmentMetadata &that) const
    {
        ISAAC_ASSERT_MSG(cluster == that.cluster && readIndex == that.readIndex,
                         "Comparison makes sense only for metadata representing the same fragment " <<
                         *this << " vs " << that);
        return
            position == that.position &&
            contigId == that.contigId &&
            reverse == that.reverse &&
            rStrandPos ==that.rStrandPos &&
            getCigarLength() == that.getCigarLength() &&
            std::equal(cigarBegin(), cigarEnd(), that.cigarBegin());
    }

    bool operator != (const FragmentMetadata &that) const
    {
        return !(*this == that);
    }

    friend std::ostream &operator<<(std::ostream &os, const FragmentMetadata &f);

    /**
     * \param readLength is needed to preallocate buffer to avoid memory operations during the processing.
     *        auxiliary code, such as unit tests, does not need to supply it.
     */
    FragmentMetadata():
        cluster(0), contigId(reference::ReferencePosition::MAX_CONTIG_ID), position(0), lowClipped(0), highClipped(0), adapterClipped_(0),
        rStrandPos(reference::ReferencePosition::NoMatch), readIndex(0), reverse(false), splitAlignment(false), decoyAlignment(false), repeatCount(0), cigarOffset(0),
        cigarLength(0), cigarBuffer(0), mismatchCount(0), uncheckedSeeds(0), gapCount(0), editDistance(0), logProbability(0.0),
        firstAnchor_(0, 0, false),
        lastAnchor_(0, 0, false),
        alignmentScore(-1U),
        mapQ(UNKNOWN_MAPQ),
        smithWatermanScore(0)
    {
    }

    FragmentMetadata(const Cluster *cluster, unsigned readIndex):
        FragmentMetadata()
    {
        this->cluster = cluster;
        this->readIndex = readIndex;
        this->lastAnchor_ = Anchor(getReadLength(), getReadLength(), false);
    }

    FragmentMetadata(
        const Cluster *cluster, unsigned readIndex,
        const unsigned short lowClipped, const unsigned short highClipped,
        const bool reverse,
        const unsigned contigId, const int64_t position,
        const bool decoyAlignment):
            FragmentMetadata(cluster, readIndex)
    {
        this->contigId = contigId;
        this->position = position;
        this->decoyAlignment = decoyAlignment;
        this->lowClipped = lowClipped;
        this->highClipped = highClipped;
        this->reverse = reverse;
    }

    FragmentMetadata(
        const Cluster *cluster, unsigned readIndex,
        const unsigned short lowClipped, const unsigned short highClipped,
        const bool reverse,
        const unsigned contigId, const int64_t position,
        const bool decoyAlignment,
        const unsigned uncheckedSeeds):
            FragmentMetadata(cluster, readIndex, lowClipped, highClipped, reverse, contigId, position, decoyAlignment)
    {
        this->uncheckedSeeds = uncheckedSeeds;
    }

    bool isReverse() const {return reverse;}
    unsigned getReadLength() const
    {
        assert(0 != cluster);
        assert(cluster->size() > readIndex);
        return (*cluster)[readIndex].getLength();
    }
    unsigned getReadIndex() const {return readIndex;}
    unsigned getObservedLength() const
    {
        ISAAC_ASSERT_MSG(!splitAlignment, "Split alignments can't have observed length" << *this);
        return rStrandPos.getPosition() - position;
    }
    unsigned getAlignmentScore() const {return alignmentScore;}
    void setAlignmentScore(unsigned as) {alignmentScore = as;}
    unsigned getCigarLength() const {return cigarLength;}
    /// Position of the first base of the fragment
    reference::ReferencePosition getFStrandReferencePosition() const
    {
        return !isNoMatch() ?
            reference::ReferencePosition(contigId, position) :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    /// Position of the last base of the fragment
    reference::ReferencePosition getRStrandReferencePosition() const
    {
        return rStrandPos - 1;
//        // ensure that the position is positive!
//        return !isNoMatch() ?
//            // observedLength can be zero if the CIGAR is soft-clipped to death or in case of
//            // split alignment with alignment position immediately following the alignment position of the last base
//            // think local translocations (most of them are fake though)
//            reference::ReferencePosition(contigId, position + std::max(observedLength, 1U) - 1) :
//            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }
    /// Position of the first aligned cycle
    reference::ReferencePosition getStrandReferencePosition() const {
        return isReverse() ? getRStrandReferencePosition() : getFStrandReferencePosition();
    }

    /// Same as f-strand position
    reference::ReferencePosition getBeginReferencePosition() const
    {
        return getFStrandReferencePosition();
    }

    /// Different from r-strand position in that it always points to the base following the last unclipped base of the fragment
    reference::ReferencePosition getEndReferencePosition() const
    {
        return !isNoMatch() ?
            rStrandPos :
            reference::ReferencePosition(reference::ReferencePosition::NoMatch);
    }

/// First cycle of fragment bcl data
    BclClusters::const_iterator getBclData() const {
        return cluster->getBclData(getReadIndex());
    }
    /// Cluster of the fragment
    const Cluster &getCluster() const {
        return *cluster;
    }

    const Read &getRead() const {
        return getCluster()[getReadIndex()];
    }

    Cigar::const_iterator cigarBegin() const
    {
        ISAAC_ASSERT_MSG(isAligned(), "Requesting CIGAR of unaligned fragment is not allowed " << *this);
        return cigarBuffer->begin() + cigarOffset;
    }

    Cigar::const_iterator cigarEnd() const
    {
        ISAAC_ASSERT_MSG(isAligned(), "Requesting CIGAR of unaligned fragment is not allowed " << *this);
        return cigarBuffer->begin() + cigarOffset + cigarLength;
    }

    uint32_t getBeginClippedLength() const
    {
        if (cigarBuffer && cigarLength)
        {
            Cigar::Component operation = Cigar::decode(cigarBuffer->at(cigarOffset));
            if (Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
        }
        return 0;
    }

    uint32_t getEndClippedLength() const
    {
        if (cigarBuffer && cigarLength)
        {
            Cigar::Component operation = Cigar::decode(cigarBuffer->at(cigarOffset + cigarLength - 1));
            if (Cigar::SOFT_CLIP == operation.second)
            {
                return operation.first;
            }
        }
        return 0;
    }

    bool isPerfectMatch() const
    {
        return !getMismatchCount() && !getBeginClippedLength() && !getEndClippedLength();
    }

    /// Unlike the observed length, excludes gaps (deletions and insertion bases)
    unsigned getMappedLength() const
    {
        ISAAC_ASSERT_MSG(cigarBuffer && cigarLength, "Read must have a valid CIGAR");
        return Cigar::getMappedLength(cigarBuffer->begin() + cigarOffset,
                                      cigarBuffer->begin() + cigarOffset + cigarLength);
    }

    /**
     * \brief Returns unadjusted position if it is adjusted due to a soft clipping
     */
    int64_t getUnclippedPosition() const
    {
        return position - getBeginClippedLength();
    }

    unsigned getMismatchCount() const {
        return mismatchCount;
    }

    bool possiblySemialigned(const reference::ContigList &contigList, const unsigned anchorLength) const
    {
        ISAAC_ASSERT_MSG(!splitAlignment, "Not checking split alignments for being semialigned:" << *this);
        const Sequence & strandSequence = getStrandSequence();
        const reference::Contig &contig = contigList.at(getContigId());

        const unsigned semialignedMismatchesMin = anchorLength / 3;

        // assume semialigned is a bad gapped aligned or an ungapped aligned with 30% mismatches on one of the ends

        return
            (gapCount && semialignedMismatchesMin < mismatchCount) ||
            makeForwardAnchor<true>(
                contig, anchorLength, semialignedMismatchesMin,
                getPosition(), strandSequence.begin(),
                getBeginClippedLength(), strandSequence.end() - getEndClippedLength()).empty() ||
            makeReverseAnchor<true>(
                contig, anchorLength, semialignedMismatchesMin, rStrandPos.getPosition(),
                strandSequence, getEndClippedLength(), strandSequence.rend() - getBeginClippedLength()).empty();
    }

    unsigned getGapCount() const {
        return gapCount;
    }

    unsigned getEditDistance() const {
        return editDistance;
    }

    ConstMismatchCycleIterator getMismatchCyclesBegin() const {return ConstMismatchCycleIterator(mismatchCycles, 1);}
    ConstMismatchCycleIterator getMismatchCyclesEnd() const {return ConstMismatchCycleIterator(mismatchCycles, MAX_CYCLES + 1);}

    void addMismatchCycle(const unsigned cycle)
    {
        ISAAC_ASSERT_MSG(cycle > 0, "Cycle numbers expected to be 1-based." << *this);
        ISAAC_ASSERT_MSG(MAX_CYCLES >= cycle, "Cycle number is too high. Check MAX_CYCLES." << *this);
        mismatchCycles.set(cycle - 1);
    }

    std::string getCigarString() const
    {
        if (cigarBuffer && cigarLength)
        {
            return Cigar::toString(*cigarBuffer, cigarOffset, cigarLength);
        }
        else
        {
            return "UNALIGNED";
        }
    }

    std::ostream &serializeCigar(std::ostream &os) const
    {
        if (cigarBuffer && cigarLength)
        {
            return Cigar::toStream(cigarBuffer->begin() + cigarOffset, cigarBuffer->begin() + cigarOffset + cigarLength, os);
        }
        else
        {
            return os << "UNALIGNED";
        }
    }

    /**
     ** \brief The cigarLength can be used to identify if a fragment has been
     ** successfully aligned
     **/
    bool isAligned() const {return 0 != cigarLength;}
    /**
     ** \brief Marks read as unaligned.
     **/
    void setUnaligned() {cigarBuffer = 0; cigarLength = 0; alignmentScore = UNKNOWN_ALIGNMENT_SCORE; mapQ = UNKNOWN_MAPQ; reverse = false, splitAlignment = false;}

    /**
     ** \brief Marks read as something that has no match position. This is different from setUnaligned
     **        as unaligned shadows still have a position of their orphan
     **/
    void setNoMatch()
        {setUnaligned(); contigId = reference::ReferencePosition::MAX_CONTIG_ID; position = 0;}
    bool isNoMatch() const {return reference::ReferencePosition::MAX_CONTIG_ID == contigId;}

    /**
     * \brief Notion of uniquely aligned in CASAVA means that a fragment was
     *        seen to have only one candidate alignment position. As this is
     *        highly dependent on the choice of seeds, alignment score based
     *        approximation should do.
     */
    bool isUniquelyAligned() const {return isAligned() && hasAlignmentScore() && (1 == repeatCount || isUnique(getAlignmentScore()));}

    bool isRepeat() const {return isAligned() && hasAlignmentScore() && REPEAT_ALIGNMENT_SCORE >= getAlignmentScore();}

    unsigned getQuality() const
    {
        const std::vector<char> &quality = getRead().getForwardQuality();

        return std::accumulate(quality.begin(), quality.end(), 0U);
    }

    typedef std::vector<char> Sequence;
    const Sequence &getStrandSequence() const {return getRead().getStrandSequence(reverse);}
    const Sequence &getStrandQuality() const {return getRead().getStrandQuality(reverse);}

    bool hasAlignmentScore() const {return UNKNOWN_ALIGNMENT_SCORE != alignmentScore;}
    bool hasMapQ() const {
        // can't assert as dodgyMapq depends on the command line argument.
        //ISAAC_ASSERT_MSG(UNKNOWN_MAPQ == mapQ || isAligned(), "mapq in unaligned: " << *this);
        return UNKNOWN_MAPQ != mapQ;}

    void incrementClipLeft(const unsigned short bases) {position += bases; if (reverse) {highClipped += bases;} else {lowClipped += bases;}}
    void incrementClipRight(const unsigned short bases) {if (!rStrandPos.isNoMatch()){ rStrandPos -= bases; } if (reverse) {lowClipped += bases;} else {highClipped += bases;}}
    void incrementAdapterClip(const unsigned short bases) {adapterClipped_ += bases;}

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short leftClipped() const {return reverse ? highClipped : lowClipped;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short rightClipped() const {return reverse ? lowClipped : highClipped;}

    /// \return number of bases clipped on the left side of the fragment in respect to the reference
    unsigned short &leftClipped() {return reverse ? highClipped : lowClipped;}
    /// \return number of bases clipped on the right side of the fragment in respect to the reference
    unsigned short &rightClipped() {return reverse ? lowClipped : highClipped;}

    void resetAlignment()
    {
        *this = FragmentMetadata(cluster, readIndex, lowClipped, highClipped, reverse, contigId, getUnclippedPosition(), decoyAlignment, uncheckedSeeds);
    }

    void resetClipping()
    {
        ISAAC_ASSERT_MSG(!isAligned(), "Alignment must be reset before clipping");
        lowClipped = 0;
        highClipped = 0;
        adapterClipped_ = 0;
    }

    unsigned getContigId() const {return contigId;}

    int64_t getPosition() const {return position;}

    unsigned updateAlignment(
        const bool collectMismatchCycles,
        const AlignmentCfg &cfg,
        const flowcell::ReadMetadata &readMetadata,
        const reference::ContigList &contigList,
        bool reverse,
        unsigned contigId,
        const int64_t strandPosition,
        const Cigar &cigarBuffer,
        const unsigned cigarOffset,
        const unsigned cigarLength = 0);

    static double calculateLogProbability(
        unsigned length,
        reference::Contig::const_iterator currentReference,
        std::vector<char>::const_iterator currentSequence,
        std::vector<char>::const_iterator currentQuality);
private:

    double calculateInsertionLogProbability(
        unsigned length,
        std::vector<char>::const_iterator currentQuality) const;

    void addMismatchCycles(
        std::vector<char>::const_iterator currentSequence,
        reference::Contig::const_iterator currentReference,
        unsigned sequenceOffset,
        unsigned length, bool reverse, const unsigned lastCycle,
        const unsigned firstCycle);

    int64_t processBackDels(const AlignmentCfg &cfg, const Cigar &cigarBuffer, const unsigned cigarEnd, unsigned &i);
    unsigned processAlign(
        const AlignmentCfg &cfg,
        unsigned &currentBase,
        const unsigned length,
        const std::vector<char>::const_iterator sequenceBegin,
        const std::vector<char>::const_iterator qualityBegin,
        const reference::Contig::const_iterator referenceBegin,
        const bool collectMismatchCycles,
        const flowcell::ReadMetadata &readMetadata,
        const bool currentReverse,
        int64_t &currentPosition);

    void processInsertion(
        const AlignmentCfg& cfg,
        unsigned &offset,
        const unsigned bases,
        std::vector<char>::const_iterator qualityBegin);

    void processDeletion(
        const AlignmentCfg& cfg,
        unsigned offset,
        const unsigned length,
        const std::vector<char>::const_iterator qualityBegin,
        const std::vector<char>::const_iterator qualityEnd,
        int64_t &currentPosition);

    void processNegativeDeletion(
        const AlignmentCfg& cfg,
        unsigned offset,
        const unsigned length,
        const std::vector<char>::const_iterator qualityBegin,
        const std::vector<char>::const_iterator qualityEnd,
        int64_t &currentPosition);

    void processSoftClip(
        unsigned &offset,
        const unsigned length,
        const std::vector<char>::const_iterator qualityBegin);

    void processHardClip(
        unsigned &offset,
        const unsigned length);

    void processContigChange(
        const AlignmentCfg& cfg,
        const unsigned newContigId,
        unsigned &currentCigarOffset,
        const reference::ContigList& contigList,
        reference::Contig::const_iterator& referenceBegin,
        unsigned &currentContigId,
        int64_t &currentPosition);

    void processFlip(
        const AlignmentCfg& cfg,
        unsigned & currentBase,
        const unsigned length,
        bool &currentReverse,
        const Read& read,
        std::vector<char>::const_iterator& sequenceBegin,
        std::vector<char>::const_iterator& qualityBegin,
        std::vector<char>::const_iterator& qualityEnd);

    template <bool firstOnly>
    Anchor makeForwardAnchor(
        const reference::Contig &contig,
        const unsigned anchorLength,
        const unsigned anchorMismatches,
        const int64_t startPosition,
        const Sequence::const_iterator &forwardBegin,
        uint32_t startOffset,
        const Sequence::const_iterator &forwardEnd) const
    {
        Sequence::const_iterator s = forwardBegin + startOffset;
        if (std::distance(s, forwardEnd) >= anchorLength)
        {
            reference::Contig::const_iterator r = contig.begin() + startPosition;
            std::size_t mismatches = countMismatches(s, r, contig.end(), anchorLength, [](char c){return c;});

            if (!firstOnly)
            {
                while (mismatches > anchorMismatches && forwardEnd != s + anchorLength && contig.end() != r + anchorLength)
                {
//                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(getCluster().getId(), "fmismatches:" << mismatches << " at " <<
//                        std::distance(forwardBegin, s) << " s:" << *s << " r:" << *r);
                    mismatches += isMismatch(*(s + anchorLength), *(r + anchorLength));
                    mismatches -= isMismatch(*s, *r);
                    ++s;
                    ++r;
                }
                startOffset = std::distance(forwardBegin, s);
            }

            if (anchorMismatches >= mismatches)
            {
                const Anchor ret(startOffset, startOffset + anchorLength, false);
//                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(getCluster().getId(), "returning fmismatches:" << mismatches << " at " <<
//                    std::distance(forwardBegin, s) << " s:" << *s << " r:" << *r << " " << ret << " " << *this);

                return ret;
            }
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(getCluster().getId(), "returning empty fanchor:" << mismatches << " at " <<
//                std::distance(forwardBegin, s) << " s:" << *s << " r:" << *r << " startPosition:" << startPosition);
        }
        // return empty anchor
        return Anchor(startOffset, startOffset, false);
    }

    template <bool firstOnly>
    Anchor makeReverseAnchor(
        const reference::Contig &contig,
        const unsigned anchorLength,
        const unsigned anchorMismatches,
        const int64_t endPosition,
        const Sequence &reverseSequence,
        uint32_t endOffset,
        const Sequence::const_reverse_iterator &reverseEnd) const
    {
        Sequence::const_reverse_iterator s = reverseSequence.rbegin() + endOffset;
        if (std::distance(s, reverseEnd) >= anchorLength)
        {
            reference::Contig::const_reverse_iterator r(contig.begin() + endPosition);
            std::size_t mismatches = countMismatches(s, r, contig.rend(), anchorLength, [](char c){return c;});

            if (!firstOnly)
            {
                while (mismatches > anchorMismatches && reverseEnd != s + anchorLength && contig.rend() != r + anchorLength)
                {
//                    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(getCluster().getId(), "rmismatches:" << mismatches << " at " <<
//                        std::distance(reverseSequence.rbegin(), s) << " s:" << *s << " r:" << *r);
                    mismatches += isMismatch(*(s + anchorLength), *(r + anchorLength));
                    mismatches -= isMismatch(*s, *r);
                    ++s;
                    ++r;
                }
                endOffset = std::distance(reverseSequence.rbegin(), s);
            }

            if (anchorMismatches >= mismatches)
            {
                const Anchor ret(reverseSequence.size() - endOffset - anchorLength, reverseSequence.size() - endOffset, false);
//                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(getCluster().getId(), "returning rmismatches:" << mismatches << " at " <<
//                    std::distance(reverseSequence.rbegin(), s) << " s:" << *s << " r:" << *r << " " << ret << " " << *this);
                return ret;
            }
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(getCluster().getId(), "returning empty ranchor:" << mismatches << " at " <<
//                std::distance(reverseSequence.rbegin(), s) << " s:" << *s << " r:" << *r);
        }
        // return empty anchor
        return Anchor(reverseSequence.size() - endOffset, reverseSequence.size() - endOffset, false);
    }

    /**
     * \brief @head - whether to update the anchor at the start or end of the sequence (in the direction of the cycles)
     */
    template <bool head,  bool firstOnly, unsigned anchorMismatches>
    bool recomputeAnchor(const reference::ContigList& contigList, const unsigned anchorLength)
    {
        if (head)
        {
            headAnchor() = computeAnchor<true, firstOnly, anchorMismatches>(contigList, anchorLength);
            return !headAnchor().empty();
        }
        else
        {
            tailAnchor() = computeAnchor<false, firstOnly, anchorMismatches>(contigList, anchorLength);
            return !tailAnchor().empty();
        }
    }


};

#ifndef _GLIBCXX_DEBUG
BOOST_STATIC_ASSERT(sizeof(FragmentMetadata) <= 256);
#endif

typedef std::vector<FragmentMetadata, common::NumaAllocator<FragmentMetadata, common::numa::defaultNodeLocal> > FragmentMetadataList;
typedef FragmentMetadataList::const_iterator FragmentIterator;

template <bool gapped>
bool putBestOnTop(
    FragmentMetadataList &fragments)
{
    FragmentMetadata &top = fragments.front();
    FragmentMetadataList::iterator best = std::min_element(fragments.begin(), fragments.end(), gapped ? FragmentMetadata::bestGappedLess : FragmentMetadata::bestUngappedLess);
    ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(top.getCluster().getId(), "    putBest" << (gapped ? "Gapped" : "Ungapped") << "OnTop\nbest: " << *best << "\ntop:" << top);
    if (fragments.begin() != best)
    {
        std::swap(*best, top);
        return true;
    }
    return false;
}


inline bool clippedByReference(const reference::ContigList& contigList, FragmentMetadata const &fragment)
{
    ISAAC_ASSERT_MSG(fragment.rStrandPos.getPosition() <= contigList.at(fragment.rStrandPos.getContigId()).size(),
                     "rStrandPos outside contig end " << fragment);
    return // ignore clipped by reference begin
        (fragment.getBeginClippedLength() && 0 == fragment.position) ||
        // ingnore clipped by reference end
        (fragment.getEndClippedLength() && fragment.rStrandPos.getPosition() == contigList.at(fragment.rStrandPos.getContigId()).size());
}


inline std::ostream &operator<<(std::ostream &os, const FragmentMetadata &f)
{
    os << "FragmentMetadata(";
    if (f.cluster)
    {
        os << common::makeFastIoString(f.getCluster().nameBegin(), f.getCluster().nameEnd()) << ",";
    }
    os << (f.cluster ? f.getCluster().getId() : 0UL) << "id "
              << f.contigId << ":"
              << f.position << ", "
              << f.rStrandPos << "rsp "
              << f.readIndex
              << (f.reverse ? 'R' : 'F') << " "
              << f.mismatchCount << "mm "
              << f.gapCount << "g "
              << f.editDistance << "ed ";
    return f.serializeCigar(os) << " "
              << f.decoyAlignment << "dcy "
              << f.splitAlignment << "sa "
              << f.logProbability << "lp "
              << f.alignmentScore << "sm "
              << int(f.mapQ) << "mq "
              << f.smithWatermanScore << "sws "
              << f.firstAnchor_<< "fa "
              << f.lastAnchor_<< "la "
              << f.leftClipped() << "lc "
              << f.rightClipped() << "rc "
              << f.uncheckedSeeds << "rs "
              << f.repeatCount << "rec)";
}


} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_METADATA_HH
