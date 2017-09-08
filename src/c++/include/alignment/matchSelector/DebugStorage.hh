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
 ** \file DebugFragmentStorage.hh
 **
 ** \brief Compares alignment result with truth data in the read name and produces various statistics for accuracy debugging.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_FRAGMENT_STORAGE_HH
#define iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_FRAGMENT_STORAGE_HH

#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/matchSelector/FragmentStorage.hh"
#include "alignment/matchSelector/debugStorage/MapqStatistics.hh"
#include "alignment/matchSelector/debugStorage/QqStatistics.hh"
#include "alignment/RestOfGenomeCorrection.hh"
#include "demultiplexing/BarcodePathMap.hh"

namespace isaac
{
namespace alignment
{


enum PairClass
{
    unaligned,
    normal,
    anomalous,
    singleton,
    maxClass = singleton,
    classesCount = maxClass + 1,
};

static const char *pairClassNames[] = {"unaligned", "normalous", "anomalous", "singleton"};

enum ReadFlags
{
    ok = 0,
    lowSm = 1,
    lowUsm = 2,
    repeatSeeds = 4,
    decoy = 8,
    flagMask = ok|lowSm|lowUsm|repeatSeeds|decoy,
    flagBits = 4,
    flagsCount = flagMask + 1
};

struct PairClassCounts
{
    std::atomic<uint64_t> totalR1Candidates_;
    std::atomic<uint64_t> decoyR1Candidates_;
    std::array<std::atomic<uint64_t>, flagsCount * flagsCount * classesCount + 1> counts_;

    ~PairClassCounts()
    {
        const uint64_t total = std::accumulate(&counts_[0], &counts_[counts_.size() - 1], uint64_t(0));
        for (std::size_t i = 0; i < counts_.size(); ++i)
        {
            const uint64_t count = counts_[i];
            if (count)
            {
//                ISAAC_THREAD_CERR << pairClassNames[getPairClass(i)] << "|" << readFlagNames[getReadAFlags(i)] << "|" << readFlagNames[getReadBFlags(i)] << ": " << count << std::endl;
                ReadFlags readAFlags = getReadAFlags(i);
                ReadFlags readBFlags = getReadBFlags(i);
                std::cerr <<
                    std::hex << std::setw(2) << i << std::dec << " " <<
                    pairClassNames[getPairClass(i)] << "|" <<
                    ((readBFlags & decoy) ? "d|" : " |") <<
                    ((readAFlags & decoy) ? "d|" : " |") <<
                    ((readBFlags & repeatSeeds) ? "rs|" : "  |") <<
                    ((readBFlags & lowUsm) ? "usm|" : "   |") <<
                    ((readBFlags & lowSm) ? "sm|" : "  |") <<
                    ((readAFlags & repeatSeeds) ? "rs|" : "  |") <<
                    ((readAFlags & lowUsm) ? "usm|" : "   |") <<
                    ((readAFlags & lowSm) ? "sm|" : "  |") <<
                    ": " << (count * 100 / total) << "% " << count << std::endl;
            }
        }

        std::cerr << "totalR1Candidates_:" << totalR1Candidates_ << std::endl;
        std::cerr << "decoyR1Candidates_:" << decoyR1Candidates_ << std::endl;
    }

    unsigned countIn(
        const BamTemplate &bamTemplate,
        const RestOfGenomeCorrection &rog,
        const bool unalgnd)
    {
        const unsigned ret = getIndex(bamTemplate, rog, unalgnd);
        ++counts_[ret];
        return ret;
    }

    static unsigned getIndex(
        const BamTemplate &bamTemplate,
        const RestOfGenomeCorrection &rog,
        const bool unalgnd)
    {
        PairClass pairClass = unalgnd ? unaligned :
            bamTemplate.isSingletonShadow() ? singleton :
            bamTemplate.isProperPair() ? normal : anomalous;

        ReadFlags r1Flags = unalgnd ? ok : getReadFlags(bamTemplate, 0, rog);
        ReadFlags r2Flags = unalgnd ? ok : getReadFlags(bamTemplate, 1, rog);

        const ReadFlags rAFlags = std::min(r1Flags, r2Flags);
        const ReadFlags rBFlags = std::max(r1Flags, r2Flags);

        const unsigned index = getIndex(rAFlags, rBFlags, pairClass);
        return index;
    }

private:
    static ReadFlags getReadFlags(
        const BamTemplate &bamTemplate,
        const unsigned readIndex,
        const RestOfGenomeCorrection &rog)
    {
        const FragmentMetadata &fragment = bamTemplate.getFragmentMetadata(readIndex);
        ReadFlags ret = ok;
        if (fragment.isAligned())
        {
            ret =  ReadFlags(ret | (fragment.alignmentScore <= 3 ? lowSm : ok));
            const unsigned usm = computeAlignmentScore(rog.getReadRogCorrection(readIndex), exp(fragment.logProbability), 0.0);
            ret = ReadFlags(ret | (usm <= 3 ? lowUsm : ok));
            ret = ReadFlags(ret | (2 < fragment.uncheckedSeeds ? repeatSeeds : ok));
            ret = ReadFlags(ret | (fragment.decoyAlignment ? decoy : ok));
        }

        return ret;
    }

    static unsigned getIndex(const ReadFlags rAFlags, const ReadFlags rBFlags, const PairClass pairClass)
    {
        return unsigned(rAFlags) | (unsigned(rBFlags) << flagBits) | ((unsigned(pairClass) << flagBits) << flagBits);
    }

    static PairClass getPairClass(std::size_t index)
    {
        return PairClass((index >> flagBits) >> flagBits);
    }

    static ReadFlags getReadAFlags(std::size_t index)
    {
        return ReadFlags(index & flagMask);
    }

    static ReadFlags getReadBFlags(std::size_t index)
    {
        return ReadFlags((index >> flagBits) & flagMask);
    }
};


namespace matchSelector
{

namespace bfs = boost::filesystem;

class DebugStorage: public FragmentStorage
{
    static const unsigned READS_MAX = 2;
    const reference::ContigList &contigList_;
    const flowcell::Layout &flowcell_;
    const AlignmentCfg &alignmentCfg_;
    const boost::filesystem::path outputDirectory_;
    const demultiplexing::BarcodePathMap barcodeQQR1Paths_;
    const demultiplexing::BarcodePathMap barcodeQQR2Paths_;
    const demultiplexing::BarcodePathMap barcodeQQR1OriPaths_;
    const demultiplexing::BarcodePathMap barcodeQQR2OriPaths_;
    const demultiplexing::BarcodePathMap barcodeMapqPaths_;
    const demultiplexing::BarcodePathMap barcodeSMPaths_;
    const demultiplexing::BarcodePathMap barcodeASPaths_;
    boost::ptr_vector<debugStorage::MapqStatistics> pairMapqStatistics_;
    boost::ptr_vector<debugStorage::QqStatistics> r1QqOriginalStatistics_;
    boost::ptr_vector<debugStorage::QqStatistics> r2QqOriginalStatistics_;
    boost::ptr_vector<debugStorage::QqStatistics> r1QqStatistics_;
    boost::ptr_vector<debugStorage::QqStatistics> r2QqStatistics_;
    alignment::matchSelector::FragmentStorage &actualStorage_;
    std::vector<std::vector<Cigar> > originalCigars_;
public:
    DebugStorage(
        const reference::ContigList &contigList,
        const AlignmentCfg &alignmentCfg,
        const flowcell::FlowcellLayoutList &flowcellLayoutList,
        const boost::filesystem::path &outputDirectory,
        const flowcell::BarcodeMetadataList &barcodeMetadataList,
        const unsigned threads,
        alignment::matchSelector::FragmentStorage &actualStorage);

    ~DebugStorage();

    virtual void store(
        const BamTemplate &bamTemplate,
        const unsigned barcodeIdx,
        const unsigned threadNumber);

    virtual void reset(const uint64_t clusterId, const bool paired)
    {
        actualStorage_.reset(clusterId, paired);
    }

    virtual void prepareFlush() noexcept
    {
        actualStorage_.prepareFlush();
    }
    virtual void flush()
    {
        actualStorage_.flush();
    }
    virtual void resize(const uint64_t clusters)
    {
        actualStorage_.resize(clusters);
    }
    virtual void reserve(const uint64_t clusters)
    {
        actualStorage_.reserve(clusters);
    }
    virtual void close()
    {
        actualStorage_.close();
    }

private:
    void updateMapqStats(
            const BamTemplate& bamTemplate,
            const unsigned barcodeIdx);
    bool restoreOriginal(
        const unsigned threadNumber,
        const std::size_t readNumber,
        FragmentMetadata &fragment);
};

} // namespace matchSelector
} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_MATCH_SELECTOR_DEBUG_FRAGMENT_STORAGE_HH
