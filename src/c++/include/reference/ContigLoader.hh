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
 ** \file ContigsLoader.hh
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_REFERENCE_CONTIGS_LOADER_HH
#define iSAAC_REFERENCE_CONTIGS_LOADER_HH

#include <boost/format.hpp>

#include "common/Threads.hpp"
#include "reference/Contig.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

void loadContig(
    const reference::SortedReferenceMetadata::Contig &xmlContig,
    ContigList::UpdateRange &contig);

template <typename ShouldLoadF> void loadContigsParallel(
    ShouldLoadF &shouldLoad,
    std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator &nextContigToLoad,
    const std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator contigsEnd,
    reference::ContigList &contigList,
    boost::mutex &mutex)
{
    const unsigned traceStep = pow(10, int(log10((contigList.size() + 99) / 100)));
    boost::lock_guard<boost::mutex> lock(mutex);
    while (contigsEnd != nextContigToLoad)
    {
        const std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator ourContig = nextContigToLoad++;
        if (shouldLoad(*ourContig))
        {
            common::unlock_guard<boost::mutex> unlock(mutex);
            const reference::SortedReferenceMetadata::Contig &xmlContig = *ourContig;
            ContigList::UpdateRange rwContig = contigList.getUpdateRange(ourContig->index_);
            loadContig(xmlContig, rwContig);
            if (!(xmlContig.index_ % traceStep))
            {
                ISAAC_THREAD_CERR << (boost::format("Contig(%3d:%8d) %s : %s\n") % xmlContig.index_ % xmlContig.totalBases_ % xmlContig.name_ % xmlContig.filePath_).str();
//                static const std::vector<char>::size_type maxBasesToPrintFromEachEnd(35);
//                ISAAC_THREAD_CERR
//                << std::string(rwContig.begin(), rwContig.begin() + std::min(rwContig.getLength()/2, maxBasesToPrintFromEachEnd)) +
//                        (rwContig.getLength() <= maxBasesToPrintFromEachEnd * 2 ? "" : " ... ") +
//                        std::string(rwContig.end() - std::min(rwContig.getLength()/2, maxBasesToPrintFromEachEnd), rwContig.end())
//                    << std::endl;
            }
        }
    }
}

/**
 * \brief loads the fasta file contigs into memory on multiple threads unless shouldLoad(contig->index_) returns false
 */
template <typename ShouldLoadF> reference::ContigList loadContigs(
    const reference::SortedReferenceMetadata::Contigs &xmlContigs,
    const std::size_t spacing,
    ShouldLoadF shouldLoad,
    common::ThreadVector &loadThreads)
{
    reference::ContigList ret(xmlContigs, spacing);
    std::vector<reference::SortedReferenceMetadata::Contig>::const_iterator nextContigToLoad = xmlContigs.begin();
    boost::mutex mutex;
    loadThreads.execute(boost::bind(&loadContigsParallel<ShouldLoadF>,
                                    boost::ref(shouldLoad),
                                    boost::ref(nextContigToLoad),
                                    xmlContigs.end(),
                                    boost::ref(ret),
                                    boost::ref(mutex)));

//    ISAAC_TRACE_STAT("loadContigs(xmlContigs) done ");

    return ret;
}

/**
 * \brief loads the fasta file contigs into memory on multiple threads
 */
template <typename AllowLoadContigT, typename IsDecoyT> reference::ContigLists loadContigs(
    const reference::SortedReferenceMetadataList &sortedReferenceMetadataList,
    const std::size_t spacing,
    const AllowLoadContigT &allowLoadContig,
    const IsDecoyT &isDecoy,
    common::ThreadVector &&loadThreads)
{
    ISAAC_TRACE_STAT("loadContigs ");

    reference::ContigLists ret;
    ret.reserve(sortedReferenceMetadataList.size());

    for(const reference::SortedReferenceMetadata &sortedReferenceMetadata : sortedReferenceMetadataList)
    {
        SortedReferenceMetadata::Contigs decoysMarkedContigs = sortedReferenceMetadata.getContigs();
        std::for_each(decoysMarkedContigs.begin(), decoysMarkedContigs.end(),
                      [&isDecoy](SortedReferenceMetadata::Contig &contig){contig.decoy_ = isDecoy(contig.name_);});

        ContigList contigList = loadContigs(decoysMarkedContigs, spacing, allowLoadContig, loadThreads);

        const std::size_t decoys =
            std::count_if(contigList.begin(), contigList.end(), [](const ContigList::Contig &contig){return contig.isDecoy();});
        ISAAC_THREAD_CERR << "Loaded " << contigList.size() << " contigs of which " << decoys << " are decoys" << std::endl;
//        ISAAC_TRACE_STAT("loadContigs(sortedReferenceMetadataList) loaded contigList ");
        ret.push_back(std::move(contigList));
//        ISAAC_TRACE_STAT("loadContigs(sortedReferenceMetadataList) pushed contigList ");
    }

    ISAAC_TRACE_STAT("loadContigs(sortedReferenceMetadataList) done ");

    return ret;
}


} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_CONTIGS_LOADER_HH
