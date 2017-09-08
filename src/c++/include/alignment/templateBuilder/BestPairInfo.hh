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
 ** \file BestPairInfo.hh
 **
 ** \brief Internal structures used by TemplateBuilder
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_BEST_PAIR_INFO_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_BEST_PAIR_INFO_HH

#include <boost/function_output_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "../../common/StaticVector.hh"
#include "common/Debug.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

static const unsigned TRACKED_REPEATS_MAX_ONE_READ = 1000;

/// arrays temporary used in buildDisjointTemplate and rescueShadow
class ShadowProbability
{
    reference::ReferencePosition pos_;
    reference::ReferencePosition rpos_;
    double logProbability_;
public:
    ShadowProbability() : logProbability_(0.0)
    {
    }
    ShadowProbability(const FragmentMetadata &shadow) :
        pos_(shadow.getFStrandReferencePosition()),
        rpos_(shadow.getRStrandReferencePosition()),
        logProbability_(shadow.logProbability)
    {
        pos_.setReverse(shadow.isReverse());
    }

    reference::ReferencePosition pos() const {return pos_;}
    double logProbability() const {return logProbability_;}
    double probability() const {return exp(logProbability_);}

    bool operator < (const ShadowProbability &that) const
    {
        // make sure same alignments stay together so that we can use std::unique to remove duplicates
        // if we happen to have fragment and its inversion, let's have the higher probability on top so that it counts towards the total
        return
            logProbability_ < that.logProbability_ ||
            (logProbability_ == that.logProbability_ &&
                (pos_ < that.pos_ ||
                    (pos_ == that.pos_ &&
                        (rpos_ < that.rpos_)
                    )
                )
            );
    }

    bool operator == (const ShadowProbability &that) const
    {
        return pos_ == that.pos_ && logProbability_ == that.logProbability_ &&
            rpos_ == that.rpos_;
    }

    bool operator != (const ShadowProbability &that) const { return !(*this == that); }

    friend std::ostream & operator << (std::ostream &os, const ShadowProbability &sp)
    {
        return os << "sp(" << sp.pos_ << ":" << sp.logProbability_<< ":" << sp.rpos_ << ")";
    }

};

struct PairProbability
{
    ShadowProbability r1_;
    ShadowProbability r2_;

    PairProbability(){}
    PairProbability(const FragmentMetadata &r1,
                    const FragmentMetadata &r2) : r1_(r1), r2_(r2){}

    bool operator < (const PairProbability &that) const
    {
        if (logProbability() < that.logProbability())
        {
            return true;
        }
        if (logProbability() != that.logProbability())
        {
            return false;
        }
        // keep same alignments together to be able to remove duplicates
        // TODO: not particularly efficient as it compares individual probabilities which we have already deal with as a sum
        return r1_ < that.r1_ || (r1_ == that.r1_ && r2_ < that.r2_);
    }

    bool operator == (const PairProbability &that) const
    {
        return r1_ == that.r1_ && r2_ == that.r2_;
    }

    bool operator != (const PairProbability &that) const {return !(*this == that);}

    double logProbability() const {return r1_.logProbability() + r2_.logProbability();}
    double probability() const {return exp(logProbability());}

    friend std::ostream & operator << (std::ostream &os, const PairProbability &pp)
    {
        return os << "PairProbability(" << pp.r1_ << "," << pp.r2_<< ")";
    }
};

struct PairInfo
{
    PairInfo() : logProbability_(-std::numeric_limits<double>::max()), swScore_(-1), matchModel_(false) {}
    PairInfo(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead, bool matchModel):
        logProbability_(oneRead.logProbability + anotherRead.logProbability),
        swScore_(oneRead.smithWatermanScore + anotherRead.smithWatermanScore), matchModel_(matchModel) {}

    void clear()
    {
        *this = PairInfo();
    }

    double probability() const {return exp(logProbability_);}

    double logProbability_;
    uint64_t swScore_;
    bool matchModel_;

    bool isWorseThan(
        const unsigned minScore,
        const RestOfGenomeCorrection &rog,
        const PairInfo &that) const
    {
        if (matchModel_== that.matchModel_)
        {
            return ISAAC_LP_LESS(logProbability_, that.logProbability_) ||
                (ISAAC_LP_EQUALS(logProbability_, that.logProbability_) && (that.swScore_ < swScore_));
        }

        const unsigned ourScoreGivenAlt = computeAlignmentScore(rog.getRogCorrection(), exp(logProbability_), exp(that.logProbability_));
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(0,"PairInfo::isWorseThan ourScoreGivenAlt: " << ourScoreGivenAlt << " minScore:" << minScore);
        if (minScore <= ourScoreGivenAlt)
        {
            return false;
        }

        const unsigned altScoreGivenUs = computeAlignmentScore(rog.getRogCorrection(), exp(that.logProbability_), exp(logProbability_));
//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(0,"PairInfo::isWorseThan altScoreGivenUs: " << altScoreGivenUs << " minScore:" << minScore);
        if (minScore <= altScoreGivenUs)
        {
            return true;
        }

//        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(0,"PairInfo::isWorseThan returning : " << !matchModel_);

        return !matchModel_;
    }

    inline bool isAsGood(const PairInfo &that) const
    {
        return that.swScore_ == swScore_ &&
            ISAAC_LP_EQUALS(that.logProbability_, logProbability_) &&
            (matchModel_ == that.matchModel_);
    }

    friend std::ostream & operator << (std::ostream & os, const PairInfo& info)
    {
        return os << "PairInfo(" << info.logProbability_ << "lp " << info.swScore_ << "sws)";
    }
};

template <typename RecordType, typename ContainerType>
void pushProbability(const RecordType &probability, ContainerType &container)
{
    if (container.end() != std::find(container.begin(), container.end(), probability))
    {
        return;
    }
    container.push_back(probability);
    // keep worst on top of the heap
    std::push_heap(container.begin(), container.end(), [](const RecordType &left, const RecordType &right){return right < left;});
    if (container.size() == container.capacity())
    {
        std::pop_heap(container.begin(), container.end(), [](const RecordType &left, const RecordType &right){return right < left;});
        // delete worst one
        container.pop_back();
    }
}

struct BestPairInfo
{
    typedef FragmentMetadataList::const_iterator FragmentIterator;

    BestPairInfo(
        const unsigned repeatThreshold = 0,
        const unsigned maxSeedsPerRead = 0)
    {
        /// Max number of orphans times the max number of shadows that can be rescued for each orphan plus the max number of shadows that have been discovered during seed matching
        reserve(repeatThreshold, maxSeedsPerRead);
    }

    void reserve(const unsigned repeatThreshold, const unsigned maxRescuedShadows)
    {
        repeats_.reserve(TRACKED_REPEATS_MAX_ONE_READ * READS_IN_A_PAIR);
    }

    void clear()
    {
        info_.clear();

        repeats_.clear();

        readProbabilities_[0].clear();
        readProbabilities_[1].clear();

        properPairProbabilities_.clear();
        allPairProbabilities_.clear();
    }

    bool empty() const
    {
        return repeats_.empty();
    }

    void appendBest(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead, const bool matchModel, const bool properPair)
    {
        ISAAC_ASSERT_MSG(oneRead.getReadIndex() != anotherRead.getReadIndex(), "Read indices match in the pair:\n" << oneRead << "\n" << anotherRead);
        if (oneRead.getReadIndex())
        {
            // make sure first is always r1
            appendBest(anotherRead, oneRead, matchModel, properPair);
            return;
        }
        if (repeats_.capacity() != repeats_.size())
        {
            const BamTemplate bamTemplate(oneRead, anotherRead, properPair);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(bamTemplate.getCluster().getId(), "appendBest:\n" << bamTemplate);
            repeats_.push_back(bamTemplate);
        }

        if (matchModel)
        {
            pushProbability(templateBuilder::PairProbability(oneRead, anotherRead), properPairProbabilities_);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oneRead.getCluster().getId(), "appendBestProperPairProbability " << templateBuilder::PairProbability(oneRead, anotherRead));
        }
        pushProbability(templateBuilder::PairProbability(oneRead, anotherRead), allPairProbabilities_);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oneRead.getCluster().getId(), "appendBestPairProbability " << templateBuilder::PairProbability(oneRead, anotherRead));
    }

    void appendPairProbability(const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead, bool matchModel)
    {
        if (oneRead.getReadIndex())
        {
            // make sure first is always r1
            appendPairProbability(anotherRead, oneRead, matchModel);
            return;
        }
        ISAAC_ASSERT_MSG(oneRead.getReadIndex() != anotherRead.getReadIndex(), "Read indices match in the pair:\n" << oneRead << "\n" << anotherRead);

        if (matchModel)
        {
            pushProbability(templateBuilder::PairProbability(oneRead, anotherRead), properPairProbabilities_);
            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oneRead.getCluster().getId(), "appendProperPairProbability " << templateBuilder::PairProbability(oneRead, anotherRead));
        }
        pushProbability(templateBuilder::PairProbability(oneRead, anotherRead), allPairProbabilities_);
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oneRead.getCluster().getId(), "appendPairProbability " << templateBuilder::PairProbability(oneRead, anotherRead));

    }

    void appendSingleProbability(const FragmentMetadata &read)
    {
        pushProbability(templateBuilder::ShadowProbability(read), readProbabilities_[read.getReadIndex()]);
    }

    void resetBest(const PairInfo &info, const FragmentMetadata &oneRead, const FragmentMetadata &anotherRead, const bool matchModel, const bool properPair)
    {
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(oneRead.getCluster().getId(), "resetBest");
        repeats_.clear();

        info_ = info;
        appendBest(oneRead, anotherRead, matchModel, properPair);
    }

    bool isWorseThan(
        const unsigned anomalousPairScoreMin,
        const RestOfGenomeCorrection &rog,
        const PairInfo &that) const
    {
        return empty() || info_.isWorseThan(anomalousPairScoreMin, rog, that);
    }

    inline bool isAsGood(const PairInfo &that) const
    {
        return !empty() && info_.isAsGood(that);
    }

    double sumUniqueReadProbabilities(const std::size_t readIndex, const double ignoreLp)
    {
        bool ignoreIgnored = false;
        double ret = 0.0;
        std::sort(readProbabilities_[readIndex].begin(), readProbabilities_[readIndex].end());

        readProbabilities_[readIndex].erase(std::unique(readProbabilities_[readIndex].begin(), readProbabilities_[readIndex].end()), readProbabilities_[readIndex].end());
        BOOST_FOREACH(const ShadowProbability &sp, readProbabilities_[readIndex])
        {
            // individual probabilities are copies. simple comparison should work
            if (!ignoreIgnored && sp.logProbability() == ignoreLp)
            {
                ignoreIgnored = true;
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(readIndex).getCluster().getId(),
                                                       "sumUniqueReadProbabilities(r" << readIndex <<
                                                       "): " << sp << " Ignored");
            }
            else
            {
                ret += sp.probability();
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(readIndex).getCluster().getId(),
                                                       "sumUniqueReadProbabilities(r" << readIndex <<
                                                       "): " << sp << " total:" << ret);
            }
        }

        return ret;
    }

    double sumUniquePairProbabilities(
        const double bestProbability, const unsigned repeatCount, const bool properOnly)
    {
        double ret = 0.0;
        ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(), "sumUniquePairProbabilities " << (!properOnly ? "all" : "proper"));

        auto &pairProbabilities = !properOnly ? allPairProbabilities_ : properPairProbabilities_;
        std::sort(pairProbabilities.begin(), pairProbabilities.end());
        pairProbabilities.erase(std::unique(pairProbabilities.begin(), pairProbabilities.end()), pairProbabilities.end());

        BOOST_FOREACH(const PairProbability &pp, pairProbabilities)
        {
            // individual probabilities are copies. simple comparison should work
            if (pp.logProbability() == bestProbability)
            {
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(),
                                                       "sumUniquePairProbabilities: " << pp << " Ignored");
            }
            else
            {
                ret += pp.probability();
                ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(),
                                                       "sumUniquePairProbabilities: " << pp << " total: " << ret);
            }
        }

        if (1 < repeatCount)
        {
            ret += exp(bestProbability) * double(repeatCount - 1);
        }

        return ret;
    }

    /**
     * \return number of unique pair repeat alignments
     */
    std::size_t removeRepeatDuplicates()
    {
        ISAAC_ASSERT_MSG(!empty(), "Invalid method call on an empty BestPairInfo" << *this);
        if (repeats_.size() == 1)
        {
            return 1;
        }

        std::sort(repeats_.begin(), repeats_.end());

        std::vector<BamTemplate>::iterator repeatsBegin = std::unique(repeats_.begin(), repeats_.end());
//        BOOST_FOREACH(const BamTemplate &repeat, std::make_pair(repeats_.begin(), repeatsBegin))
//        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(),
//                                                   "removeRepeatDuplicates keeping: " << repeat);
//        }
//
//        BOOST_FOREACH(const BamTemplate &repeat, std::make_pair(repeatsBegin, repeats_.end()))
//        {
//            ISAAC_THREAD_CERR_DEV_TRACE_CLUSTER_ID(repeats_.front().getFragmentMetadata(0).getCluster().getId(),
//                                                   "removeRepeatDuplicates removing: " << repeat);
//        }

        repeats_.erase(repeatsBegin, repeats_.end());
        return repeats_.size();
    }

    unsigned bestPairMismatchCount() const
    {
        return repeats_.front().getFragmentMetadata(0).getMismatchCount() + repeats_.front().getFragmentMetadata(1).getMismatchCount();
    }

//    bool isKUnique() const
//    {
//        return !repeats_.empty() && repeats_.front().isKUnique();
//    }

    std::size_t repeatsCount() const {return repeats_.size();}

    double probability() const {return info_.probability();}
    double logProbability() const {return info_.logProbability_;}

    const BamTemplate &repeat(std::size_t index) const {return repeats_.at(index);}
private:
    // corresponding entries represent read pairs. All read pairs in repeatAlignments_ represent alignments of the same bestTemplateScore_ to the repeat
    static const unsigned READS_IN_A_PAIR = 2;
    static const std::size_t MAX_PROBABILITIES_TO_KEEP = 3;
    PairInfo info_;
    common::StaticVector<templateBuilder::ShadowProbability, MAX_PROBABILITIES_TO_KEEP + 1> readProbabilities_[READS_IN_A_PAIR];
    common::StaticVector<templateBuilder::PairProbability, MAX_PROBABILITIES_TO_KEEP + 1> properPairProbabilities_;
    common::StaticVector<templateBuilder::PairProbability, MAX_PROBABILITIES_TO_KEEP + 1> allPairProbabilities_;

    std::vector<BamTemplate> repeats_;

    friend std::ostream & operator << (std::ostream & os, const BestPairInfo& bestPairInfo)
    {
        os << "BestPairInfo(";
        if (bestPairInfo.empty())
        {
            return os << "empty)";
        }
        if (!bestPairInfo.repeats_.empty())
        {
            os << bestPairInfo.repeats_.front() << ",";
        }
        return os << bestPairInfo.info_ << " " << bestPairInfo.bestPairMismatchCount() << "bem)";
    }

};

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_BEST_PAIR_INFO_HH
