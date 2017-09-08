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
 ** \file BinMetadata.hh
 **
 ** \brief Metadata associated to the unsorted alignment results
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_ALIGNMENT_BIN_METADATA_HH
#define iSAAC_ALIGNMENT_BIN_METADATA_HH

#include <numeric>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "flowcell/BarcodeMetadata.hh"
#include "reference/ReferencePosition.hh"
#include "io/Fragment.hh"

namespace isaac
{
namespace alignment
{

struct BinMetadata;
inline std::ostream &operator<<(std::ostream &os, const BinMetadata &binMetadata);

struct BarcodeCounts
{
    BarcodeCounts() : elements_(0), gaps_(0), splits_(0), cigarLength_(0){}
    // total number of elements in the bin barcode
    uint64_t elements_;
    // total number of gaps in the bin barcode reads
    uint64_t gaps_;
    // total number of splits in the bin barcode reads. Splits are count normally once except for FLIPs which could
    // count twice unless FLIP is not followed by CONTIG or position adjustment
    uint64_t splits_;
    // sum of all fragment cigar lengths in the bin barcode.
    uint64_t cigarLength_;

    BarcodeCounts& operator += (const BarcodeCounts &that)
    {
        elements_ += that.elements_;
        gaps_ += that.gaps_;
        splits_ += that.splits_;
        cigarLength_ += that.cigarLength_;
        return *this;
    }

    friend BarcodeCounts operator + (const BarcodeCounts &left, const BarcodeCounts &right)
    {
        BarcodeCounts ret(left);
        return ret += right;
    }
};

class BinMetadata
{
    unsigned binIndex_;
    /// first genomic position covered by the bin
    reference::ReferencePosition binStart_;
    /// bin length in bases
    uint64_t length_;
    boost::filesystem::path binFilePath_;
    // offset from the beginning of the data file.
    // Note that single file can later be broken down into multiple BinMetadata objects
    uint64_t dataOffset_;
    // number of bytes stored in binFilePath_ at dataOffset_
    uint64_t dataSize_;
    uint64_t seIdxElements_;
    uint64_t rIdxElements_;
    uint64_t fIdxElements_;
    uint64_t nmElements_;
    std::vector<BarcodeCounts> barcodeBreakdown_;

    /*
     * \brief enable serialization
     */
    template <class Archive> friend void serialize(Archive &ar, BinMetadata &bm, const unsigned int version);

public:
    BinMetadata() :
        binIndex_(0),
        binStart_(0),
        length_(0),
        dataOffset_(0),
        dataSize_(0),
        seIdxElements_(0),
        rIdxElements_(0),
        fIdxElements_(0),
        nmElements_(0),
        barcodeBreakdown_(){}

    BinMetadata(
        const unsigned barcodesCount,
        const unsigned binIndex,
        const reference::ReferencePosition binStart,
        const uint64_t length,
        const boost::filesystem::path &binFilepath) :
            binIndex_(binIndex),
            binStart_(binStart),
            length_(length),
            binFilePath_(binFilepath),
            dataOffset_(0),
            dataSize_(0),
            seIdxElements_(0),
            rIdxElements_(0),
            fIdxElements_(0),
            nmElements_(0),
            barcodeBreakdown_(barcodesCount){}

    void swap(BinMetadata &that) throw()
    {
        using std::swap;
        swap(binIndex_, that.binIndex_);
        swap(binStart_, that.binStart_);
        swap(length_, that.length_);
        swap(binFilePath_, that.binFilePath_);
        swap(dataOffset_, that.dataOffset_);
        swap(dataSize_, that.dataSize_);
        swap(seIdxElements_, that.seIdxElements_);
        swap(rIdxElements_, that.rIdxElements_);
        swap(fIdxElements_, that.fIdxElements_);
        swap(nmElements_, that.nmElements_);
        swap(barcodeBreakdown_, that.barcodeBreakdown_);
    }

    friend void swap(BinMetadata &left, BinMetadata &right) throw()
    {
        left.swap(right);
    }

    void merge(const BinMetadata &that)
    {
        if (!isUnalignedBin())
        {
            ISAAC_ASSERT_MSG(binStart_ + length_ <= that.binStart_, "Unexpected order of bins in merge: " << *this << " and " << that);
            length_ = that.binStart_ + that.length_ - binStart_;
        }
        else
        {
            ISAAC_ASSERT_MSG(that.isUnalignedBin() && dataOffset_ + dataSize_ == that.dataOffset_, "Attempt to merge non adjacent bins: " << *this << " and " << that);
        }

        dataSize_ += that.dataSize_;
        seIdxElements_ += that.seIdxElements_;
        rIdxElements_ += that.rIdxElements_;
        fIdxElements_ += that.fIdxElements_;
        nmElements_ += that.nmElements_;
        std::transform(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), that.barcodeBreakdown_.begin(), barcodeBreakdown_.begin(), std::plus<BarcodeCounts>());
    }

    unsigned getIndex() const
    {
        return binIndex_;
    }

    reference::ReferencePosition getBinStart() const
    {
        return binStart_;
    }

    bool hasPosition(const reference::ReferencePosition pos) const
    {
        return binStart_ <= pos && pos < (binStart_ + length_);
    }

    reference::ReferencePosition getBinEnd() const
    {
        return isUnalignedBin() ? reference::ReferencePosition(reference::ReferencePosition::NoMatch) : binStart_ + length_;
    }

    bool coversPosition(const reference::ReferencePosition pos) const
    {
        ISAAC_ASSERT_MSG(!isUnalignedBin(), "Checking positions in unaligned bins is not allowed");
        return getBinStart() <= pos && getBinEnd() > pos;
    }

    bool isUnalignedBin() const {return binStart_.isTooManyMatch();}

    const boost::filesystem::path & getPath() const
    {
        return binFilePath_;
    }
    const std::string &getPathString() const
    {
        return binFilePath_.string();
    }

    bool sameContig(const BinMetadata &that) const
    {
        // direct comparison of path objects leads to another dynamic string allocation
        return (isUnalignedBin() && that.isUnalignedBin()) ||
            (isUnalignedBin() == that.isUnalignedBin() && binStart_.getContigId() == that.binStart_.getContigId());
    }

    bool samePath(const BinMetadata &that) const
    {
        // direct comparison of path objects leads to another dynamic string allocation
        return binFilePath_.string() == that.binFilePath_.string();
    }

    uint64_t getDataOffset() const
    {
        return dataOffset_;
    }

    uint64_t getDataSize() const
    {
        return dataSize_;
    }

    bool isEmpty() const
    {
        return 0 == dataSize_;
    }

    uint64_t getLength() const {return length_;}

    /**
     * \brief increment the the corresponding chunk size and total data size.
     */
    void incrementDataSize(const reference::ReferencePosition pos, const uint64_t by)
    {
        dataSize_ += by;
    }

    /**
     * \brief increment the the corresponding chunk size and total data size.
     */
    void incrementDataSize(const uint64_t recordNumber, const uint64_t by)
    {
        dataSize_ += by;
    }

    uint64_t getSeIdxElements() const
    {
        return seIdxElements_;
    }

    void incrementSeIdxElements(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).elements_ += by;
        seIdxElements_ += by;
    }

    uint64_t getRIdxElements() const
    {
        return rIdxElements_;
    }

    void incrementRIdxElements(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).elements_ += by;
        rIdxElements_ += by;
    }

    uint64_t getFIdxElements() const
    {
        return fIdxElements_;
    }

    void incrementFIdxElements(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).elements_ += by;
        fIdxElements_ += by;
    }

    uint64_t getNmElements() const
    {
        return nmElements_;
    }

    void incrementNmElements(const uint64_t sequenceHash, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).elements_ += by;
        nmElements_ += by;
    }

    void incrementGapCount(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).gaps_ += by;
    }

    void incrementSplitCount(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).splits_ += by;
    }

    void incrementCigarLength(const reference::ReferencePosition pos, const uint64_t by, const unsigned barcodeIdx)
    {
        barcodeBreakdown_.at(barcodeIdx).cigarLength_ += by;
    }

    uint64_t getTotalElements() const
    {
        //return  getSeIdxElements() + getRIdxElements() + getFIdxElements() + getNmElements();
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::elements_, _2)));
    }

    uint64_t getBarcodeElements(const unsigned barcodeIdx) const
    {
        return barcodeBreakdown_.at(barcodeIdx).elements_;
    }

    uint64_t getBarcodeGapCount(const unsigned barcodeIdx) const
    {
        ISAAC_ASSERT_MSG(barcodeBreakdown_.size() > barcodeIdx, "Invalid barcode requested: " << barcodeIdx << " for size: " << barcodeBreakdown_.size());
        return barcodeBreakdown_[barcodeIdx].gaps_;
    }

    uint64_t getTotalGapCount() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::gaps_, _2)));
    }

    uint64_t getTotalSplitCount() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::splits_, _2)));
    }

    uint64_t getTotalCigarLength() const
    {
        return std::accumulate(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), 0,
                               boost::bind(std::plus<uint64_t>(),
                                           _1, boost::bind(&BarcodeCounts::cigarLength_, _2)));
    }

    uint64_t getDataEndOffset() const
    {
        return getDataOffset() + getDataSize();
    }

    /**
     * \brief restarts counting from current offset
     */
    void startNew()
    {
        dataOffset_ += dataSize_;
        dataSize_ = 0;
        seIdxElements_ = 0;
        rIdxElements_ = 0;
        fIdxElements_ = 0;
        nmElements_ = 0;
        std::fill(barcodeBreakdown_.begin(), barcodeBreakdown_.end(), BarcodeCounts());
    }
};


struct BinMetadataList : public std::vector<alignment::BinMetadata>
{
    BinMetadataList(){}
    BinMetadataList(size_t size) : std::vector<alignment::BinMetadata>(size){}
    BinMetadataList(size_t size, const BinMetadata &bm) : std::vector<alignment::BinMetadata>(size, bm){}
};
typedef boost::reference_wrapper<BinMetadata> BinMetadataCRef;
typedef std::vector<BinMetadataCRef >BinMetadataCRefList;


inline std::ostream &operator<<(std::ostream &os, const BinMetadata &binMetadata)
{
    return os << "BinMetadata("
              << binMetadata.getIndex() << "id "
              << binMetadata.getBinStart() << "bs "
              << binMetadata.getLength() << "bl "
              << binMetadata.getDataSize() << "ds "
              << binMetadata.getDataOffset() << "do "
              << binMetadata.getSeIdxElements() << "se "
              << binMetadata.getRIdxElements() << "rs "
              << binMetadata.getFIdxElements() << "f "
              << binMetadata.getPathString() << ")";
}

} // namespace alignment
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_BIN_METADATA_HH
