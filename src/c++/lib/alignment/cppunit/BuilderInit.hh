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

#include <vector>
#include <string>
#include <boost/foreach.hpp>

#include "alignment/BclClusters.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"
#include "reference/Contig.hh"

inline isaac::flowcell::ReadMetadataList getReadMetadataList()
{
    const unsigned l0 = 100, l1 = 100;
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, l0, 0, 0))
            (isaac::flowcell::ReadMetadata(l0 + 1, l0 + l1, 1, l0))
            ;
    return ret;
}

inline isaac::flowcell::ReadMetadataList getReadMetadataList(const unsigned l0, const unsigned l1)
{
    std::vector<isaac::flowcell::ReadMetadata> ret =
        boost::assign::list_of
            (isaac::flowcell::ReadMetadata(1, l0, 0, 0))
            (isaac::flowcell::ReadMetadata(l0 + 1, l0 + l1, 1, l0))
            ;
    return ret;
}

inline std::string getContig(const std::string name, const unsigned length)
{
    char bases[] = {'A', 'C', 'G', 'T'};
    std::string contig;
    contig.resize(length);
    BOOST_FOREACH(char &base, contig)
    {
        base = bases[rand() % 4];
    }
    return contig;
}

inline void show(const std::vector<char> &s)
{
    BOOST_FOREACH(char c, s) std::cerr << c;
}

template <typename ContainerT>
inline std::vector<char> reverseComplement(const ContainerT &forward)
{
    std::vector<char> tr(256, 'N');
    tr['A'] = 'T';
    tr['C'] = 'G';
    tr['G'] = 'C';
    tr['T'] = 'A';
    std::vector<char> reverse;
    reverse.reserve(forward.size());
    for (typename ContainerT::const_reverse_iterator b = forward.rbegin(); forward.rend() != b; ++b)
    {
        reverse.push_back(tr[*b]);
    }
    return reverse;
}

inline std::vector<char> vectorFromString(const std::string &str)
{
    return std::vector<char>(str.begin(), str.end());
}

inline std::vector<char> vectorFromString(const std::string &str, bool ignoreSpaces)
{
    if (ignoreSpaces)
    {
        std::vector<char> ret;
        std::remove_copy(str.begin(), str.end(), std::back_inserter(ret), ' ');
        return ret;
    }
    return std::vector<char>(str.begin(), str.end());
}

template <typename ContainerT>
inline std::string substr(const ContainerT &from, std::string::size_type __pos = 0,
                   std::string::size_type __n = std::string::npos)
{
    return std::string(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

template <typename ContainerT>
inline std::vector<char> subv(const ContainerT &from, std::string::size_type __pos = 0,
                       std::string::size_type __n = std::string::npos)
{
    return std::vector<char>(from.begin() + __pos, std::string::npos == __n ? from.end() : from.begin() + __pos + __n);
}

//inline std::vector<char> subv(const std::vector<char> &from, std::string::size_type __pos,
//                       std::string::size_type __n)
//{
//    return std::vector<char>(from.begin() + __pos, from.begin() + __pos + __n);
//}

//template <typename ContainerT>
//inline std::vector<char> subv(const ContainerT &from, std::string::size_type __pos)
//{
//    return std::vector<char>(from.begin() + __pos, from.end());
//}

inline std::vector<char> operator +(const std::vector<char> &right, const isaac::reference::Contig &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}


inline std::vector<char> operator +(const std::vector<char> &right, const std::vector<char> &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}

inline std::vector<char> operator +(const std::vector<char> &right, const std::string &left)
{
    std::vector<char> ret(right);
    ret.insert(ret.end(), left.begin(), left.end());
    return ret;
}


struct TestContigList : public isaac::reference::ContigList
{
    TestContigList(TestContigList &&that) : isaac::reference::ContigList(std::move(that)){}
    TestContigList(){}

    TestContigList(const std::vector<char> &reference) : TestContigList(std::vector<std::vector<char> >(1, reference))
    {
    }

    TestContigList(const std::string &reference) : TestContigList(std::vector<std::string>(1, reference))
    {
    }

    template <typename RefT>
    static isaac::reference::SortedReferenceMetadata makeSortedReferenceMetadata(const std::vector<RefT> &contigs)
    {
        isaac::reference::SortedReferenceMetadata sortedReferenceMetadata;
        std::size_t genomicOffset = 0;
        for (const RefT &contig: contigs)
        {
            sortedReferenceMetadata.putContig(
                genomicOffset, "chr" + std::to_string(sortedReferenceMetadata.getContigsCount() + 1), "blah.fa",
                genomicOffset, contig.size(),
                contig.size(),
                contig.size(), sortedReferenceMetadata.getContigsCount(),
                "", "", "");
            genomicOffset += contig.size();
        }
        return sortedReferenceMetadata;
    }

    template <typename RefT>
    TestContigList(const std::vector<RefT> &contigs) : isaac::reference::ContigList(makeSortedReferenceMetadata(contigs).getContigs(), 1000)
    {
        for (std::size_t contigId = 0; contigId < contigs.size(); ++contigId)
        {
            isaac::reference::ContigList::UpdateRange rwContig = getUpdateRange(contigId);
            std::copy(contigs[contigId].begin(), contigs[contigId].end(), rwContig.begin());
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

template <typename ContainerT>
std::vector<char> getBcl(const ContainerT &bases)
{
    std::vector<char> bcl;
    bcl.reserve(bases.size());
    BOOST_FOREACH(char b, bases)
    {
        using isaac::oligo::getValue;
        bcl.push_back((40 << 2) | getValue(b));
    }
    return bcl;
}

inline std::vector<char> getBclVector(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const isaac::reference::ContigList &contigList,
    const unsigned contigId, const int offset0, const int offset1,
    const bool reverse0 = false, const bool reverse1 = true)
{
    const isaac::reference::Contig &contig = contigList[contigId];
    const unsigned length0 = readMetadataList[0].getLength();
    const unsigned length1 = readMetadataList[1].getLength();

    std::string forward;
    std::string reverse;
    forward.reserve(contig.size());
    reverse.reserve(contig.size());
    BOOST_FOREACH(const char base, contig)
    {
        using isaac::oligo::getReverseBase;
        using isaac::oligo::getValue;
        forward.push_back(base);
        reverse.push_back(getReverseBase(getValue(base)));
    }

    std::reverse(reverse.begin(), reverse.end());
    const std::string &s0 = reverse0 ? reverse : forward;
    const std::string &s1 = reverse1 ? reverse : forward;
    std::string bases(s0.begin() + offset0, s0.begin() + offset0 + length0);
    bases += std::string(s1.begin() + offset1, s1.begin() + offset1 + length1);
    return getBcl(bases);
}

inline isaac::alignment::BclClusters getBclClusters(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const std::vector<char> &clusters)
{
    isaac::alignment::BclClusters ret(isaac::flowcell::getTotalReadLength(readMetadataList));
    ret.reserveClusters(1, false);
    std::copy(clusters.begin(), clusters.end(), ret.cluster(0));
    return ret;
}

inline isaac::alignment::BclClusters getBcl(
    const std::vector<isaac::flowcell::ReadMetadata> &readMetadataList,
    const isaac::reference::ContigList &contigList,
    const unsigned contigId, const int offset0, const int offset1,
    const bool reverse0 = false, const bool reverse1 = true)
{
    const std::vector<char> clusters = getBclVector(
        readMetadataList, contigList, contigId, offset0, offset1, reverse0, reverse1);
    return getBclClusters(readMetadataList, clusters);
}
