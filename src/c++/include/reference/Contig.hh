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
 ** \file Contig.hh
 **
 ** \brief Definition of a contig
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_CONTIG_HH
#define iSAAC_REFERENCE_CONTIG_HH

#include <iostream>
#include <string>
#include <vector>

#include "common/Debug.hh"
#include "common/NumaContainer.hh"
#include "common/SameAllocatorVector.hh"
#include "reference/ReferencePosition.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{
template <typename AllocatorT>
struct BasicReferenceSequence : public std::vector<char, AllocatorT>
{
    typedef std::vector<char, AllocatorT> BaseT;
    BasicReferenceSequence(const AllocatorT &a) : BaseT(a){}
    BasicReferenceSequence() : BaseT(){}
    explicit BasicReferenceSequence(
        typename BaseT::size_type n, const typename BaseT::value_type& value = typename BaseT::value_type(),
        const typename BaseT::allocator_type& a = typename BaseT::allocator_type()) :
        BaseT(n, value, a){}

    template <typename OtherT>
    BasicReferenceSequence &operator =(const OtherT &that) {this->assign(that.begin(), that.end()); return *this;}
};

template <typename AllocatorT>
struct BasicContig
{
public:
    typedef BasicReferenceSequence<AllocatorT> ReferenceSequence;
    typedef typename ReferenceSequence::const_iterator const_iterator;
    typedef typename ReferenceSequence::const_reverse_iterator const_reverse_iterator;

    const char &back() const {return *(end_ - 1);}
    const_iterator begin() const {return begin_;}
    bool empty() const {return !size();}
    const_iterator end() const {return end_;}

    char operator[](const std::size_t offset) const {return *(begin_ + offset);}
    void swap(BasicContig &that)
    {
        std::swap(index_, that.index_);
        std::swap(name_, that.name_);
        std::swap(decoy_, that.decoy_);
        std::swap(begin_, that.begin_);
        std::swap(end_, that.end_);
    }

    friend void swap(BasicContig &left, BasicContig &right)
    {
        left.swap(right);
    }

    const_reverse_iterator rbegin() const {return const_reverse_iterator(end_);}
    const_reverse_iterator rend() const {return const_reverse_iterator(begin_);}

    std::size_t size() const {return std::distance(begin_, end_);}


    const std::string &getName() const {return name_;}

    BasicContig(const unsigned index, const std::string &name, const bool decoy, const_iterator begin, const_iterator end) :
        index_(index), name_(name), decoy_(decoy), begin_(begin), end_(end){;}

    bool isDecoy() const {return decoy_;}

    friend std::ostream &operator <<(std::ostream &os, const BasicContig &contig)
    {
        return os << "Contig(" << contig.index_ << "," << contig.name_ << "," << (contig.decoy_ ?  "dcy," : ",") << contig.size() << ")";
    }

    unsigned getIndex() const {return index_;}

private:
    unsigned index_;
    std::string name_;
    bool decoy_;
    const_iterator begin_;
    const_iterator end_;
};


template <typename AllocatorT>
class BasicContigList : protected std::vector<BasicContig<AllocatorT > >
{
    static const std::size_t CONTIG_LENGTH_MIN = 0x10000;
    static const std::size_t OFFSET_MAX = 0x0ffffffffUL;
    static const std::size_t TRANSLATION_TABLE_SIZE = (OFFSET_MAX + 1) / CONTIG_LENGTH_MIN;

//    typedef BasicContigList<AllocatorT> MyT;
//    typedef typename AllocatorT::template rebind<char> CharAllocatorRebind;
//    typedef typename CharAllocatorRebind::other CharAllocator;

public:

    typedef BasicContig<AllocatorT> Contig;
    typedef std::vector<Contig> BaseT;

    typedef BasicReferenceSequence<AllocatorT> ReferenceSequence;
    typedef typename ReferenceSequence::iterator ReferenceSequenceIterator;
    typedef typename ReferenceSequence::const_iterator ReferenceSequenceConstIterator;
    typedef unsigned short ContigId;
    typedef uint32_t Offset;
    static const ContigId INVALID_CONTIG_ID = ContigId(0) - 1;

protected:
    std::vector<ContigId> contigIdFromScaledOffset_;
    ReferenceSequence referenceSequence_;

public:

    using BaseT::size;
    using typename BaseT::const_iterator;
    using typename BaseT::value_type;
    typedef AllocatorT allocator_type;

    const Contig &at(std::size_t contigIndex) const {return BaseT::at(contigIndex);}
    const Contig &front() const {return BaseT::front();}
    const Contig &operator [](std::size_t contigIndex) const {return BaseT::operator[](contigIndex);}
    const_iterator begin() const {return BaseT::begin();}
    Offset beginOffset(const std::size_t contigIndex) const {return std::distance(referenceSequence_.begin(), at(contigIndex).begin());}
    const_iterator end() const {return BaseT::end();}
    Offset endOffset(const std::size_t contigIndex) const {return std::distance(referenceSequence_.begin(), at(contigIndex).end());}
    Offset endOffset() const {return std::distance(referenceSequence_.begin(), referenceSequence_.end());}

    ReferenceSequenceConstIterator referenceBegin() const {return referenceSequence_.begin();}

    struct UpdateRange : std::pair<ReferenceSequenceIterator, ReferenceSequenceIterator>
    {
        typedef std::pair<ReferenceSequenceIterator, ReferenceSequenceIterator> BaseT;
        UpdateRange(ReferenceSequenceIterator begin, ReferenceSequenceIterator end) : BaseT(begin, end){}
        ReferenceSequenceIterator begin() {return this->first;}
        ReferenceSequenceIterator end() {return this->second;}
        std::size_t getLength() const {return std::distance(this->first, this->second);}
        friend std::ostream &operator << (std::ostream &os, const UpdateRange &ur)
        {
            return os << "UpdateRange(length=" << ur.getLength() << ")";
        }
    };
    UpdateRange getUpdateRange(const std::size_t contigIndex)
    {
        const Contig &contig = at(contigIndex);
        return UpdateRange(referenceSequence_.begin() + std::distance<ReferenceSequenceConstIterator>(referenceSequence_.begin(), contig.begin()),
                           referenceSequence_.begin() + std::distance<ReferenceSequenceConstIterator>(referenceSequence_.begin(), contig.end()));
    }

    ContigId contigIdFromOffset(const Offset offset) const {return contigIdFromScaledOffset_.at(offset / CONTIG_LENGTH_MIN);}
    std::int64_t positionFromOffset(const Offset offset) const {return offset - contigBeginOffset(contigIdFromOffset(offset));}
//
//    isaac::reference::ReferencePosition referencePositionFromOffset(const Offset offset, bool reverse) const
//    {
//        const ContigId contigId = contigIdFromOffset(offset);
//        return isaac::reference::ReferencePosition(contigId, offset - contigBeginOffset(contigId), false, reverse);
//    }

    BasicContigList() : BaseT(){}
    BasicContigList(const isaac::reference::SortedReferenceMetadata::Contigs &xmlContigs, const std::size_t spacing);

    BasicContigList(const BasicContigList &that): BaseT()
    {
//        ISAAC_THREAD_CERR << "BasicContigList copy constructor" << std::endl;
        assign(that);
//        ISAAC_THREAD_CERR << "BasicContigList copy constructor done" << std::endl;
    }
    BasicContigList(const BasicContigList &that, const AllocatorT &a): BaseT(),
        contigIdFromScaledOffset_(that.contigIdFromScaledOffset_), referenceSequence_(a)
    {
        assign(that);
    }
//    template <typename OtherContigLisT>
//    BasicContigList &operator = (const OtherContigLisT &that)
//    {
//        assing(that);
//        return *this;
//    }
//
    template <typename OtherContigLisT>
    void assign(const OtherContigLisT& that)
    {
        contigIdFromScaledOffset_ = that.contigIdFromScaledOffset_;
        referenceSequence_ = that.referenceSequence_;
        for (const Contig &contig : that)
        {
            this->push_back(
                Contig(contig.getIndex(), contig.getName(), contig.isDecoy(),
                       referenceSequence_.begin() + std::distance(that.referenceSequence_.begin(), contig.begin()),
                       referenceSequence_.begin() + std::distance(that.referenceSequence_.begin(), contig.end())));
        }
    }


    BasicContigList(BasicContigList &&that)
    {
        BaseT::swap(that);
        contigIdFromScaledOffset_.swap(that.contigIdFromScaledOffset_);
        referenceSequence_.swap(that.referenceSequence_);
    }

    BasicContigList &operator =(BasicContigList &&that)
    {
        if (this != &that)
        {
            contigIdFromScaledOffset_.swap(that.contigIdFromScaledOffset_);
            referenceSequence_.swap(that.referenceSequence_);
            BaseT::swap(that);
        }
        return *this;
    }

    void swap(BasicContigList &that)
    {
        contigIdFromScaledOffset_.swap(that.contigIdFromScaledOffset_);
        referenceSequence_.swap(that.referenceSequence_);
        BaseT::swap(that);
    }

    friend void swap(BasicContigList &left, BasicContigList &right)
    {
        left.swap(right);
    }

    std::size_t contigBeginOffset(const ContigId contigId) const
    {
        return std::distance(referenceSequence_.begin(), at(contigId).begin());
    }

private:
    BaseT &base() {return *this;}
};

// keep a separate copy of linear reference on each numa node
typedef BasicContigList<common::NumaAllocator<char, 0> > ContigList;
typedef common::SameAllocatorVector<ContigList, common::NumaAllocator<char, 0> > ContigLists;
class NumaContigLists
{
    common::NumaContainerReplicas<ContigLists> replicas_;
public:

    NumaContigLists(ContigLists &&node0Lists) :replicas_(std::move(node0Lists))
    {
        ISAAC_THREAD_CERR << "NumaContigLists ContigLists constructor" << std::endl;
    }

    const ContigLists &node0Container() const {return replicas_.node0Container();}
    const ContigLists &threadNodeContainer() const {return replicas_.threadNodeContainer();}
//    operator const ContigLists &()const {return replicas_.threadNodeContainer();}
};

typedef typename ContigList::Contig Contig;



/// Total length of all the contigs of a genome

inline std::size_t genomeLength(const ContigList &contigList)
{
    return std::accumulate(
        contigList.begin(), contigList.end(),
        size_t(0), boost::bind<size_t>(std::plus<size_t>(), _1, boost::bind(&ContigList::Contig::size, _2)));
}


} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIG_HH
