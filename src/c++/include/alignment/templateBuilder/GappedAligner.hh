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
 ** \file GappedAligner.hh
 **
 ** \brief Uses banded Smith-Waterman algorithm to align fragment
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_GAPPED_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_GAPPED_ALIGNER_HH

#include "alignment/templateBuilder/AlignerBase.hh"
#include "alignment/BandedSmithWaterman.hh"

namespace isaac
{
namespace alignment
{

namespace templateBuilder
{

class GappedAligner: public AlignerBase
{
public:
    GappedAligner(
        const bool collectMismatchCycles,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const bool smartSmithWaterman,
        const unsigned smithWatermanGapSizeMax,
        const AlignmentCfg &alignmentCfg);

    bool realignBadUngappedAlignments(
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        FragmentMetadataList &fragmentList,
        FragmentSequencingAdapterClipper &adapterClipper,
        Cigar &cigarBuffer);

    /**
     ** \brief Calculate the gapped alignment of a fragment
     **/
    unsigned alignGapped(
        const flowcell::ReadMetadata &readMetadata,
        const FragmentSequencingAdapterClipper &adapterClipper,
        const reference::ContigList &contigList,
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer)
    {
        return alignGapped(smartSmithWaterman_, readMetadata, adapterClipper, contigList, fragmentMetadata, cigarBuffer);
    }

    unsigned alignGapped(
        const bool smartSmithWaterman,
        const flowcell::ReadMetadata &readMetadata,
        const FragmentSequencingAdapterClipper &adapterClipper,
        const reference::ContigList &contigList,
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer);
protected:
    static const unsigned HASH_KMER_LENGTH = 7;
    static const unsigned QUERY_LENGTH_MAX = 65536;

    const bool smartSmithWaterman_;
    const unsigned smithWatermanGapSizeMax_;
    typedef BandedSmithWaterman<16> Bsw16;
    typedef BandedSmithWaterman<32> Bsw32;
    typedef BandedSmithWaterman<64> Bsw64;
    Bsw16 bandedSmithWaterman16_;
    Bsw32 bandedSmithWaterman32_;
    Bsw64 bandedSmithWaterman64_;

    common::StaticVector<unsigned, 2> hashedQueryTile_;
    common::StaticVector<unsigned, 2> hashedQueryCluster_;
    common::StaticVector<unsigned, 2> hashedQueryReadIndex_;
    common::StaticVector<std::vector<char>::const_iterator, 2> hashedQueryBegin_;
    common::StaticVector<std::vector<char>::const_iterator, 2> hashedQueryEnd_;

    // initialize all k-mers to the magic value -1 (NOT_FOUND)
    static const std::size_t MAX_KMER = oligo::MaxKmer<HASH_KMER_LENGTH, unsigned short>::value;
    common::StaticVector<unsigned short, MAX_KMER + 1 > queryKmerOffsets_[2];
    static const unsigned short UNINITIALIZED_OFFSET_MAGIC = static_cast<unsigned short>(-1);
    static const unsigned short REPEAT_OFFSET_MAGIC = static_cast<unsigned short>(-2);
    // count of hits required to assume that part of the sequence will anchor at a position.
    // the smaller the number, the more reads will go into smith-waterman. The higher the number
    // the more likely is that some good gapped alignments will not be attempted
    static const unsigned char SUFFICIENT_NUMBER_OF_HITS = 4;

    // counts of corresponding offsets of the first base of the data from the first base of the reference
    // offset of 0 is in the middle of the trackedOffsets.
    common::StaticVector<unsigned char, QUERY_LENGTH_MAX> trackedOffsets_;


    bool makesSenseToGapAlign(
        const unsigned tile, const unsigned cluster, const unsigned read, const bool reverse,
        const std::vector<char>::const_iterator queryBegin,
        const std::vector<char>::const_iterator queryEnd,
        const reference::Contig::const_iterator databaseBegin,
        const reference::Contig::const_iterator databaseEnd);

    /**
     ** \brief Calculate the gapped alignment of a fragment and accept it if it is better
     **/
    template <typename BswT>
    bool realignOne(
        BswT &bandedSmithWaterman,
        const bool smartSmithWaterman,
        const unsigned smitWatermanGapsMax,
        const reference::ContigList& contigList,
        const flowcell::ReadMetadata& readMetadata,
        FragmentMetadata& fragmentMetadata,
        templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
        Cigar& cigarBuffer);

    template <typename BswT>
    unsigned alignGapped(
        BswT &bandedSmithWaterman,
        const bool smartSmithWaterman,
        const flowcell::ReadMetadata &readMetadata,
        const FragmentSequencingAdapterClipper &adapterClipper,
        const reference::ContigList &contigList,
        FragmentMetadata &fragmentMetadata,
        Cigar &cigarBuffer);

    template <typename BswT>
    bool realignBadUngappedAlignments(
        BswT &bandedSmithWaterman,
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        FragmentMetadataList &fragmentList,
        FragmentSequencingAdapterClipper &adapterClipper,
        Cigar &cigarBuffer);


private:
    void updateComponent(const unsigned cigarOffset, uint64_t len,
                         const Cigar::OpCode op, Cigar& cigarBuffer);
};

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_GAPPED_ALIGNER_HH
