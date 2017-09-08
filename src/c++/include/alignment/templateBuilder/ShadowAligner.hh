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
 ** \file ShadowAligner.hh
 **
 ** \brief Aligns shadows
 ** 
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_SHADOW_ALIGNER_HH
#define iSAAC_ALIGNMENT_SHADOW_ALIGNER_HH

#include <vector>
#include <boost/noncopyable.hpp>

#include "reference/Contig.hh"
#include "alignment/BandedSmithWaterman.hh"
#include "alignment/Cigar.hh"
#include "alignment/Cluster.hh"
#include "alignment/SequencingAdapter.hh"
#include "alignment/TemplateLengthStatistics.hh"
#include "alignment/templateBuilder/GappedAligner.hh"
#include "alignment/templateBuilder/UngappedAligner.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

/**
 ** \brief Utility component aligning shadows.
 **
 ** The intended use case is for the TemplateBuilder to delegate the alignment
 ** of shadow reads (or poorly aligned mates) to this specialized component.
 **
 **/
template <unsigned SHADOW_KMER_LENGTH>
class ShadowAligner: boost::noncopyable
{
public:
    /**
     ** \brief Construct a ShadowAligner for a reference genome and a given
     ** set of reads.
     **/
    ShadowAligner(
        const bool collectMismatchCycles,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const unsigned gappedMismatchesMax,
        const unsigned smitWatermanGapsMax,
        const bool smartSmithWaterman,
        const bool noSmithWaterman,
        const bool splitAlignments,
        const AlignmentCfg &alignmentCfg,
        Cigar &cigarBuffer);
    /**
     ** \brief Helper method to align the shadow of an orphan.
     **
     ** Uses the position, read index, and orientation of the orphan to infers
     ** the orientation and range of positions of the shadow (according to the
     ** templateLengthStatistics_).
     **
     ** Hashes the all k-mers of length orphanKmerLength_ from the orphan
     ** sequence into the orphaKmerPositions_.
     **
     ** Find all the possible alignment positions of the orphan by comparing the
     ** k-mers of length orphanKmerLength_ from the region of the reference
     ** where the orphan could align to the hashed k-mers from the orphan
     ** sequence.
     **
     ** Align the orphan to the candidate positions.
     **
     ** \param orphan[in] the correctly aligned orphan
     **
     ** \param shadow[out] the unaligned shadow
     **/
    bool rescueShadows(
        const reference::ContigList &contigList,
        const FragmentMetadata &orphan,
        const unsigned shadowsMax,
        FragmentMetadataList &shadowList,
        const flowcell::ReadMetadata &shadowReadMetadata,
        templateBuilder::FragmentSequencingAdapterClipper &adapterClipper,
        const TemplateLengthStatistics &templateLengthStatistics);
    const Cigar &getCigarBuffer() const {return cigarBuffer_;}
private:

    /// Length of the k-mers used to rescue shadows and misaligned reads
    static const unsigned shadowKmerCount_ = (1 << (2 * SHADOW_KMER_LENGTH));
    const unsigned gappedMismatchesMax_;
    const unsigned smitWatermanGapsMax_;
    const bool noSmithWaterman_;
    const bool splitAlignments_;
    const flowcell::FlowcellLayoutList &flowcellLayoutList_;
    const templateBuilder::UngappedAligner ungappedAligner_;
    Cigar &cigarBuffer_;
    /**
     ** \brief Cached storage for the position of the k-mers in the shadow
     **
     ** Note that this is a really fast and cheap but imperfect to rescue
     ** shadows or mis-aligned reads. The index used to access elements in the
     ** vector is made from a k-mer of length shadowKmerLength_ (the vector has
     ** 4 ^ shadowKmerLength_ positions). The values in the table at position i
     ** is the first position in the read where the k-mer was found (-1 if not
     ** found). Repeats are recorded only once in the table. This allows to
     ** identify extremely quickly if a k-mer in the reference belongs to the
     ** read. The shadowKmerLength_ should stay small enough to ensure that the
     ** table stays in the L1 cache.
     **/
//    std::vector<short> shadowKmerPositions_;
    common::StaticVector<short, shadowKmerCount_> shadowKmerPositions_;
    /// Hash all the k-mers of length shadowKmerLength_ into shadowKmerPositions_
    unsigned hashShadowKmers(const std::vector<char> &sequence);
    /**
     ** \brief Cached storage for the candidate start positions of the shadow
     **
     ** Note: all these positions are relative to the beginning of the region
     ** where the shadow is expected to be (as opposed to the position relative
     ** to the beginning of the reference).
     **/
    std::vector<int64_t> shadowCandidatePositions_;
    /// Find all candidate positions for a shadow sequence on a given reference interval
    bool findShadowCandidatePositions(
        const std::pair<int64_t, int64_t> &alignmentStartPositionRange,
        const int64_t referenceOffset,
        const reference::Contig::const_iterator referenceBegin,
        const reference::Contig::const_iterator referenceEnd,
        const std::vector<char> &shadowSequence,
        std::vector<int64_t> &shadowCandidatePositions);
    bool findShadowCandidatePositions(
        const FragmentMetadata& orphan,
        const TemplateLengthStatistics& templateLengthStatistics,
        const Read& shadowRead, const bool shadowReverse,
        const reference::Contig& contig,
        std::vector<int64_t> &shadowCandidatePositions);
    bool alignCandidates(
        const FragmentMetadata& orphan,
        const bool reverse, const reference::ContigList& contigList,
        const flowcell::ReadMetadata& readMetadata,
        templateBuilder::FragmentSequencingAdapterClipper& adapterClipper,
        FragmentMetadataList& shadowList);
};

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_SHADOW_ALIGNER_HH
