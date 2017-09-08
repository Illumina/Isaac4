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
 **/

#include <vector>
#include <string>
#include <boost/foreach.hpp>

#include "alignment/BclClusters.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"
#include "reference/Contig.hh"


inline std::string getContig(const std::string name, const unsigned length)
{
    char bases[] = {'A', 'C', 'G', 'T'};
    std::string contig;
    contig.resize(length);
    for(char &base : contig)
    {
        base = bases[rand() % 4];
    }
    return contig;
}


inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

struct TestContigList : public isaac::reference::ContigList
{
    TestContigList(TestContigList &&that) : isaac::reference::ContigList(std::move(that)){}
    TestContigList(){}

    TestContigList(const std::string &reference)
    {
        referenceSequence_ = vectorFromString(reference);
        push_back(isaac::reference::Contig(0, "tada", false, referenceSequence_.begin(), referenceSequence_.end()));
    }

    template <typename RefT>
    TestContigList(const std::vector<RefT> &contigs)
    {
        referenceSequence_.reserve(std::accumulate(contigs.begin(), contigs.end(), 0UL,
                                          [](const std::size_t &a, const RefT &b){return a + b.size();}));
        for (const RefT &contig: contigs)
        {
            const std::size_t before = referenceSequence_.size();
            referenceSequence_.insert(referenceSequence_.end(), contig.begin(), contig.end());
            push_back(isaac::reference::Contig(0, "chr" + std::to_string(size() +1), false, referenceSequence_.begin() + before, referenceSequence_.end()));
        }
    }

    TestContigList &operator = (TestContigList &&that)
    {
        isaac::reference::ContigList::operator =(std::move(that));
        return *this;
    }
};

inline TestContigList getContigList(const unsigned l0 = 210, const unsigned l1 = 220, const unsigned l4 = 60)
{
    const std::string c2 = getContig("c2", 230);
    std::string v = ("AAAAA");
    v.insert(v.end(), c2.begin(), c2.end());
    const std::string c3(v);
    const std::string c4 = getContig("c4", l4);
    const std::string c1 = getContig("c1", l1);
    const std::string c0 = getContig("c0", l0);

    return TestContigList(boost::assign::list_of(c0)(c1)(c2)(c3)(c4).convert_to_container<std::vector<std::string> >());
}

