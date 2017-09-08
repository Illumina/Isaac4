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
 ** \file PeAdapterTrimmer.hh
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_PE_ADAPTER_TRIMMER_HH
#define iSAAC_ALIGNMENT_TEMPLATE_BUILDER_PE_ADAPTER_TRIMMER_HH

//#include "alignment/BamTemplate.hh"
//#include "alignment/BandedSmithWaterman.hh"
//#include "alignment/Cigar.hh"
//#include "alignment/Cluster.hh"
#include "alignment/FragmentMetadata.hh"
//#include "alignment/Match.hh"
//#include "alignment/RestOfGenomeCorrection.hh"
//#include "alignment/templateBuilder/FragmentBuilder.hh"
//#include "alignment/templateBuilder/ShadowAligner.hh"
//#include "alignment/TemplateLengthStatistics.hh"
#include "flowcell/ReadMetadata.hh"
#include "reference/Contig.hh"
//#include "templateBuilder/BestPairInfo.hh"

namespace isaac
{
namespace alignment
{
namespace templateBuilder
{

/**
 ** \brief Utility component trimming adapters in PE alignments resulting in template shorter than read length
 **     This particular alignment will result in the bases marked as A soft-clipped on read ends:
 **     <AAAAAAASSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 **             SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSAAAAAAA>
 **/
class PeAdapterTrimmer: boost::noncopyable
{
public:
    PeAdapterTrimmer(
        const bool collectMismatchCycles,
        const bool trimPEAdapters,
        const AlignmentCfg &alignmentCfg);

    bool checkTrimPEAdapter(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        FragmentMetadata &r1Fragment,
        FragmentMetadata &r2Fragment,
        Cigar &cigarBuffer) const;

    bool trimPEAdapter(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        FragmentMetadata & forwardFragment,
        FragmentMetadata & reverseFragment,
        Cigar &cigarBuffer) const;

    FragmentMetadata trimForwardPEAdapter(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const FragmentMetadata & forwardFragment,
        const reference::ReferencePosition &adapterPosition,
        Cigar &cigarBuffer) const;

    bool trimForwardPEAdapter(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        uint32_t cycles,
        FragmentMetadata & forwardFragment,
        Cigar &cigarBuffer) const;

    FragmentMetadata trimReversePEAdapter(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        const FragmentMetadata &reverseFragment,
        const reference::ReferencePosition &adapterPosition,
        Cigar &cigarBuffer) const;

    bool trimReversePEAdapter(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadata &readMetadata,
        uint32_t cycles,
        FragmentMetadata & reverseFragment,
        Cigar &cigarBuffer) const;

    bool trimPEAdapterCycles(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const unsigned cycles,
        FragmentMetadata &fragment,
        Cigar &cigarBuffer) const;

    void trimPEAdapterCycles(
        const reference::ContigList &contigList,
        const flowcell::ReadMetadataList &readMetadataList,
        const unsigned cycles,
        FragmentMetadataList &fragments,
        Cigar &cigarBuffer) const;
private:
    const bool trimPEAdapters_;
    const bool collectMismatchCycles_;
    const AlignmentCfg alignmentCfg_;
};

} // namespace templateBuilder
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_TEMPLATE_BUILDER_PE_ADAPTER_TRIMMER_HH
