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
 ** \file Match.hh
 **
 ** \brief Abstract representation of the match of a seed to a reference position.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_HH
#define iSAAC_ALIGNMENT_MATCH_HH

#include <iostream>

#include "common/Numa.hh"
#include "oligo/Kmer.hh"
#include "reference/Contig.hh"
#include "reference/ReferencePosition.hh"


namespace isaac
{
namespace alignment
{

/**
 ** \brief a component that represents a match from a seed to a given reference location.
 **
 **/
struct Match
{
    Match(reference::ContigList::Offset contigListOffset,
          bool reverse) :
              contigListOffset_(contigListOffset) , reverse_(reverse){}

    reference::ContigList::Offset contigListOffset_;
    bool reverse_;

    bool operator == (const Match &that)
    {
        return that.contigListOffset_ == contigListOffset_ && that.reverse_ == reverse_;
    }
};

typedef std::vector<Match, common::NumaAllocator<Match, common::numa::defaultNodeLocal> > Matches;
typedef std::vector<Matches> MatchLists;

inline std::ostream &operator<<(std::ostream &os, const Match &match)
{
    return os << "Match(" << match.contigListOffset_ << (match.reverse_ ? "r" : "f") << ")";
}

} //namespace alignment
} //namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_HH
