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
 ** \file BamTemplate.hh
 **
 ** \brief DNA/RNA sequence composed of one or several Fragment(s), as defined
 ** by the SAM Format Specification
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_HH
#define iSAAC_ALIGNMENT_TEMPLATE_HH

#include <iostream>
#include <numeric>
#include <vector>

#include "alignment/FragmentMetadata.hh"
#include "alignment/RestOfGenomeCorrection.hh"

namespace isaac
{
namespace alignment
{


struct PairInfo
{
    PairInfo() : logProbability_(-std::numeric_limits<double>::max()), probability_(1.0),
        swScore_(-1), properPair_(false)//, mismatchCount_(-1)
    {}
    PairInfo(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead, bool properPair_):
        logProbability_(oneRead.logProbability + anotherRead.logProbability),
        probability_(exp(logProbability_)),
        swScore_(oneRead.smithWatermanScore + anotherRead.smithWatermanScore), properPair_(properPair_)//,
//        mismatchCount_(oneRead.getMismatchCount() + anotherRead.getMismatchCount())
    {}

    void clear()
    {
        *this = PairInfo();
    }

    double probability() const {return probability_;}

    double logProbability_;
    double probability_;
    int swScore_;
    bool properPair_;
//    unsigned mismatchCount_;

    bool isWorseThan(
        const double anomalousPairHandicap,
        const double rog,
        const PairInfo &that) const
    {
        if (properPair_== that.properPair_)
        {
            return ISAAC_LP_LESS(logProbability_, that.logProbability_) ||
                (ISAAC_LP_EQUALS(logProbability_, that.logProbability_) && (that.swScore_ < swScore_));
        }

//        if (mismatchCount_ + anomalousPairHandicap <= that.mismatchCount_)
//        {
//            return false;
//        }

        if (anomalousPairHandicap <= probability_ / (that.probability_ + rog))
        {
            return false;
        }
//
//        const unsigned ourScoreGivenAlt = computeAlignmentScore(rog, exp(logProbability_), exp(that.logProbability_));
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(0,"PairInfo::isWorseThan ourScoreGivenAlt: " << ourScoreGivenAlt << " anomalousPairHandicap:" << anomalousPairHandicap);
//        if (anomalousPairHandicap <= ourScoreGivenAlt)
//        {
//            return false;
//        }

//        if (that.mismatchCount_ + anomalousPairHandicap <= mismatchCount_)
//        {
//            return true;
//        }

        if (anomalousPairHandicap <= that.probability_ / (probability_ + rog))
        {
            return true;
        }

//        const unsigned altScoreGivenUs = computeAlignmentScore(rog, exp(that.logProbability_), exp(logProbability_));
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(0,"PairInfo::isWorseThan altScoreGivenUs: " << altScoreGivenUs << " anomalousPairHandicap:" << anomalousPairHandicap);
//        if (anomalousPairHandicap <= altScoreGivenUs)
//        {
//            return true;
//        }

//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(0,"PairInfo::isWorseThan returning : " << !properPair_);

        return !properPair_;
    }

    inline bool isAsGood(const PairInfo &that) const
    {
        return that.swScore_ == swScore_ &&
            ISAAC_LP_EQUALS(that.logProbability_, logProbability_) &&
            (properPair_ == that.properPair_);
    }

    friend std::ostream & operator << (std::ostream & os, const PairInfo& info)
    {
        return os << "PairInfo(" << info.logProbability_ << "lp " << info.swScore_ << "sws" << info.properPair_ << "pp)";
    }
};


/**
 ** \brief Container to encapsulate all the data and metadata associated to a DNA/RNA BamTemplate.
 **
 ** \sa Fragment
 ** \sa FragmentId
 ** \sa Cigar
 ** \sa FragmentMetadata
 **/
class BamTemplate
{
public:
    BamTemplate();
    BamTemplate(const FragmentMetadata &read1, const FragmentMetadata &read2,
        const bool properPair = false, const unsigned alignmentScore = -1U);
    BamTemplate(const FragmentMetadata &read1);
    BamTemplate(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster);
    /**
     ** \brief Initialization for a given cluster
     **
     ** Create the appropriate unaligned FragmentMetadata for each read in the
     ** cluster.
     **/
    void reset(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster);

    unsigned getMismatchCount() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getMismatchCount, _2)));
    }

    unsigned getQuality() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0U,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getQuality, _2)));
    }

    unsigned getEditDistance() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getEditDistance, _2)));
    }

    unsigned getTotalReadLength() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::getReadLength, _2)));
    }

    bool isUnanchored() const
    {
        return fragmentMetadataList_.end() ==
            std::find_if(fragmentMetadataList_.begin(), fragmentMetadataList_.end(),
                         boost::bind(&FragmentMetadata::getAlignmentScore, _1) != 0);
    }

    bool isUniquelyAligned() const
    {
        return hasAlignmentScore() && isUnique(getAlignmentScore());
    }

    double getLogProbability() const
    {
        return pairInfo_.logProbability_;
    }

    int getSmithWatermanScore() const
    {
        return pairInfo_.swScore_;
    }

    unsigned gapCount() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), 0,
                               boost::bind(std::plus<unsigned>(), _1,
                                           boost::bind(&FragmentMetadata::gapCount, _2)));
    }

    bool isBetterThan(const double anomalousPairHandicap, const double rog, const BamTemplate &alt) const
    {
        return alt.pairInfo_.isWorseThan(anomalousPairHandicap, rog, pairInfo_);
    }

    bool isRepeat() const
    {
        return std::accumulate(fragmentMetadataList_.begin(), fragmentMetadataList_.end(), true,
                               boost::bind(std::logical_and<bool>(), _1,
                                           boost::bind(&FragmentMetadata::isRepeat, _2)));
    }

    unsigned getFragmentCount() const {return fragmentMetadataList_.size();}
    const FragmentMetadata &getFragmentMetadata(unsigned fragmentIndex) const {return fragmentMetadataList_.at(fragmentIndex);}
    const FragmentMetadata &getMateFragmentMetadata(const FragmentMetadata &mate) const {return fragmentMetadataList_.at(getFragmentCount() - 1 - mate.getReadIndex());}
    unsigned getAlignmentScore() const {return alignmentScore_;}
    void resetAlignmentScore() {alignmentScore_ = -1U;}
    bool hasAlignmentScore() const {return -1U != alignmentScore_;}
    void setAlignmentScore(unsigned alignmentScore) {alignmentScore_ = alignmentScore;}
    bool getPassesFilter() const {return fragmentMetadataList_[0].getCluster().getPf();}
    bool isProperPair() const {return pairInfo_.properPair_;}
    bool isSingletonShadow() const {return 2 == getFragmentCount() && fragmentMetadataList_[0].isAligned() !=  fragmentMetadataList_[1].isAligned();}
    bool filterLowQualityFragments(const int mapqThreshold);
    const Cluster &getCluster() const {return fragmentMetadataList_[0].getCluster();}

    void setNoMatch()
    {
        fragmentMetadataList_[0].setNoMatch();
        if (2 == getFragmentCount())
        {
            fragmentMetadataList_[1].setNoMatch();
        }
    }

    void setDodgy(const unsigned char dodgyAlignmentScore)
    {
        fragmentMetadataList_[0].mapQ = dodgyAlignmentScore;
        if (2 == getFragmentCount())
        {
            fragmentMetadataList_[1].mapQ = dodgyAlignmentScore;
        }
    }

    void setFragmentMapq(const unsigned fragmentIndex, const unsigned char mapq)
    {
        fragmentMetadataList_[fragmentIndex].mapQ = mapq;
    }

    BclClusters::const_iterator nameBegin() const
    {
        return fragmentMetadataList_[0].getCluster().nameBegin();
    }

    BclClusters::const_iterator nameEnd() const
    {
        return fragmentMetadataList_[0].getCluster().nameEnd();
    }

    unsigned getNameLength() const {return std::distance(nameBegin(), nameEnd());}

    bool operator ==(const BamTemplate &right) const
    {
        ISAAC_ASSERT_MSG(getFragmentCount() == right.getFragmentCount(), "Incorrect template comparison");

        if (getFragmentMetadata(0) != right.getFragmentMetadata(0))
        {
            return false;
        }

        if (1 == getFragmentCount())
        {
            return true;
        }

        return getFragmentMetadata(1) == right.getFragmentMetadata(1);
    }

    bool operator <(const BamTemplate &right) const
    {
        ISAAC_ASSERT_MSG(getFragmentCount() == right.getFragmentCount(), "Incorrect template comparison");
        if (getFragmentMetadata(0) < right.getFragmentMetadata(0))
        {
            return true;
        }

        if (1 == getFragmentCount())
        {
            return false;
        }

        if (right.getFragmentMetadata(0) < getFragmentMetadata(0))
        {
            return false;
        }

        return getFragmentMetadata(1) < right.getFragmentMetadata(1);
    }
//    unsigned debugClass_ = 0;
private:
    friend std::ostream &operator<<(std::ostream &os, const BamTemplate &bamTemplate);

    common::StaticVector<FragmentMetadata, 2> fragmentMetadataList_;

    /**
     ** This depends on the all the pLogCorrect values for all the possible
     ** alignments for this template across the whole reference. It also takes into
     ** account the rest-of-genome correction. Value of -1U indicates unknown alignment score.
     **/
    unsigned alignmentScore_;
    PairInfo pairInfo_;
};

inline std::ostream &operator<<(std::ostream &os, const BamTemplate &bamTemplate)
{
    return 2 == bamTemplate.getFragmentCount() ?
        os << "BamTemplate(\n"
        << bamTemplate.fragmentMetadataList_.at(0) << "\n" << bamTemplate.fragmentMetadataList_.at(1) << "," <<
        bamTemplate.alignmentScore_ << "as," << bamTemplate.pairInfo_ << ")" :
        os << "BamTemplate("
        << bamTemplate.fragmentMetadataList_.at(0) << "," <<
        bamTemplate.alignmentScore_ << "as," << bamTemplate.pairInfo_<< ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_HH
