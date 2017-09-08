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
 ** \file SplitReadAligner.hh
 **
 ** \brief Uses seed match discrepancies to detect insertions, deletions, translocations, inversions in the reads that span single event
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH
#define iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH


#include "alignment/templateBuilder/AlignerBase.hh"
#include "alignment/TemplateLengthStatistics.hh"

namespace isaac
{
namespace alignment
{

namespace templateBuilder
{

class SplitReadAligner: public AlignerBase
{
public:
    SplitReadAligner(
        const bool collectMismatchCycles,
        const AlignmentCfg &alignmentCfg);

    void alignSimpleSv(
        Cigar &cigarBuffer,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const TemplateLengthStatistics &templateLengthStatistics,
        FragmentMetadataList &fragmentList) const;

    bool resolveConflict(
        const reference::ContigList& contigList,
        const flowcell::ReadMetadata& readMetadata,
        const bool regularIndelsOnly,
        Cigar& cigarBuffer,
        FragmentMetadataList& fragmentList,
        const FragmentMetadata &head, const FragmentMetadata &tail) const;
private:

    bool alignIndel(
        Cigar &cigarBuffer,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const bool regularIndelsOnly,
        FragmentMetadata &head,
        const FragmentMetadata &tail) const;

    bool alignTranslocation(
        Cigar &cigarBuffer,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        FragmentMetadata &head,
        const FragmentMetadata &tail) const;

    bool alignSimpleDeletion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        const FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata) const;

    bool alignLeftAnchoredInversion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        const FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata) const;

    bool alignRightAnchoredInversion(
        Cigar &cigarBuffer,
        FragmentMetadata &headFragment,
        const unsigned headSeedOffset,
        const FragmentMetadata &tailFragment,
        const unsigned tailSeedOffset,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata) const;

    bool alignSimpleInsertion(
        Cigar &cigarBuffer,
        const FragmentMetadata &headAlignment,
        const unsigned headSeedOffset,
        FragmentMetadata &tailAlignment,
        const unsigned tailSeedOffset,
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata) const;

    bool mergeDeletionAlignments(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const unsigned bestMismatches,
        const int deletionLength,
        const flowcell::ReadMetadata &readMetadata) const;

    bool mergeLeftAnchoredInversions(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const unsigned bestMismatches,
        const int deletionLength,
        const flowcell::ReadMetadata &readMetadata) const;

    bool mergeRightAnchoredInversionAlignments(
        Cigar &cigarBuffer,
        FragmentMetadata &headAlignment,
        const FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const unsigned bestMismatches,
        const int deletionLength,
        const flowcell::ReadMetadata &readMetadata) const;

    bool mergeInsertionAlignments(
        Cigar &cigarBuffer,
        const FragmentMetadata &headAlignment,
        FragmentMetadata &tailAlignment,
        const unsigned bestOffset,
        const reference::ContigList &contigList,
        const unsigned bestMismatches,
        const unsigned insertionLength,
        const flowcell::ReadMetadata &readMetadata) const;
        
    bool pickBestSplit(const bool tmp1Worked, const FragmentMetadata& tmp1, const bool tmp2Worked, const FragmentMetadata& tmp2,
                       FragmentMetadataList& fragmentList) const;
};

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_FRAGMENT_BUILDER_SIMPLE_INDEL_ALIGNER_HH
