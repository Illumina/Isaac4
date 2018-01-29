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
 ** \file BamLoader.hh
 **
 ** Component to read Bam files.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BAM_BAM_PARSER_HH
#define iSAAC_BAM_BAM_PARSER_HH

#include <boost/iterator/iterator_facade.hpp>

#include "alignment/Cigar.hh"
#include "common/Endianness.hh"
#include "common/Exceptions.hh"
#include "common/FastIo.hh"
#include "flowcell/ReadMetadata.hh"
#include "oligo/Nucleotides.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace bam
{

struct BamParserException : common::IoException
{
    BamParserException(const std::string &message) : common::IoException(EINVAL, message){}
};


inline char bamBase(const unsigned char bamSeq)
{
    static const unsigned char BAM_BASES[] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    const unsigned char bamBase = BAM_BASES[bamSeq];
    return bamBase;
}

inline unsigned char bamToBcl(const unsigned char qual, const unsigned char bamSeq)
{
    static const oligo::Translator<> translator = {};
    const unsigned char q = 0xFF == qual ? 0 : std::min<unsigned char>(qual, 0x3f);
    const unsigned char base = bamBase(bamSeq);
//            std::cerr << bamBase;
//            std::cerr << char(q + 33);
    const unsigned char baseValue = translator[base];
    return oligo::INVALID_OLIGO == baseValue ? 0 : (baseValue | (q << 2));
}

struct BamBlockHeader : boost::noncopyable
{
protected:
    int32_t block_size;
    int32_t refID;
    int32_t pos;
    uint32_t bin_mq_nl;
    uint32_t flag_nc;
    int32_t l_seq;
    int32_t next_refID;
    int32_t next_pos;
    int32_t tlen;
    char read_name[1];
public:

    static const unsigned MULTI_SEGMENT = 0x01 << 16;
    static const unsigned UNMAPPED_SEGMENT = 0x04 << 16;
    static const unsigned NEXT_UNMAPPED_SEGMENT = 0x08 << 16;
    static const unsigned REV_COMPL = 0x10 << 16;
    static const unsigned NEXT_REV_COMPL = 0x20 << 16;
    static const unsigned FIRST_SEGMENT = 0x40 << 16;
    static const unsigned LAST_SEGMENT = 0x80 << 16;
    static const unsigned VERNDOR_FAILED = 0x200 << 16;
    static const unsigned SECONDARY_ALIGNMENT = 0x100 << 16;
    static const unsigned SUPPLEMENTARY_ALIGNMENT = 0x800 << 16;

    unsigned char getReadNameLength() const {return bin_mq_nl;}
    unsigned short getCigarLength() const {return flag_nc;}
    const char *nameBegin() const {return read_name ;}
    const char *nameEnd() const {return read_name + getReadNameLength();}
    const unsigned *getCigar() const {return reinterpret_cast<const unsigned*>(read_name + getReadNameLength());}
    const unsigned char *getSeq() const {return reinterpret_cast<const unsigned char *>(getCigar() + getCigarLength());}
    const unsigned char *getQual() const {return getSeq() + (getLSeq() + 1) / 2;}
    bool isPaired() const {return flag_nc & MULTI_SEGMENT;}
    bool isReverse() const {return flag_nc & REV_COMPL;}
    bool isNextReverse() const {return flag_nc & NEXT_REV_COMPL;}
    bool isReadOne() const
    {
        // FIRST_SEGMENT is not set for single-ended data
        return !(flag_nc & LAST_SEGMENT);
    }

    bool isPf() const {return !(flag_nc & VERNDOR_FAILED);}
    bool isSupplementaryAlignment() const {return flag_nc & SUPPLEMENTARY_ALIGNMENT;}
    bool isSecondaryAlignment() const {return flag_nc & SECONDARY_ALIGNMENT;}
    bool isUnmapped() const {return flag_nc & UNMAPPED_SEGMENT;}
    bool isNextUnmapped() const {return flag_nc & NEXT_UNMAPPED_SEGMENT;}

    int getRefId() const {return common::extractLittleEndian<int>(&refID);}
    int getNextRefId() const {return common::extractLittleEndian<int>(&next_refID);}

    int getPos() const {return common::extractLittleEndian<int>(&pos);}
    int getNextPos() const {return common::extractLittleEndian<int>(&next_pos);}
    int getLSeq() const {return common::extractLittleEndian<int>(&l_seq);}

    friend std::ostream &seqToStream(std::ostream &os, const BamBlockHeader &bamBlockHeader)
    {
        const unsigned char * const pSeq = bamBlockHeader.getSeq();
        unsigned seqLen = bamBlockHeader.getLSeq();

        unsigned seqOffset = 0;
        while(seqLen--)
        {
            os << bamBase((*(pSeq + seqOffset / 2) >> (4 * ((seqOffset + 1) % 2))) & 0x0F);
            ++seqOffset;
        }

        return os;
    }

    friend std::ostream &operator << (std::ostream &os, const BamBlockHeader &bamBlockHeader)
    {
        return os << "BamBlockHeader("  << bamBlockHeader.read_name <<  "," <<
            bamBlockHeader.refID << ":" << bamBlockHeader.pos <<
            ")";
//        return seqToStream(os, bamBlockHeader) << ")";
    }

#pragma pack(push, 1)
    struct Tag
    {
        char tag[2];
        char value_type;
        union
        {
            char z;
            int8_t c;
            uint8_t C;
            int16_t s;
            uint16_t S;
            int32_t i;
            uint32_t I;

            //TODO: support binary one day.
        } value;

        std::size_t valueLength() const
        {
            switch (value_type)
            {
            case 'Z':
                return strlen(&value.z) + 1;
            case 'c':
            case 'C':
                return 1;
            case 's':
            case 'S':
                return 2;
            case 'i':
            case 'I':
                return 4;
            }
            ISAAC_ASSERT_MSG(false, "Value type " << value_type << " is not supported");
            return -1;
        }

        const Tag * end() const
        {
            return reinterpret_cast<const Tag*>(&value.z + valueLength());
        }

        bool isTag(const char *str) const
        {
            return str[0] == tag[0] && str[1] == tag[1];
        }

        const char *valueBegin() const
        {
            return &value.z;
        }

        const char *valueEnd() const
        {
            return &value.z + valueLength();
        }
        friend std::ostream &operator << (std::ostream &os, const Tag &t)
        {
            return os << "Tag(" << t.tag[0] << t.tag[1] <<
                    "," << t.value_type << "," <<
                    (t.value_type == 'Z' ? &t.value.z : "integer") << ")";
        }
    };
#pragma pack(pop)

    class TagIterator :
        public boost::iterator_facade<
        TagIterator
        , const Tag &
        , boost::forward_traversal_tag
        >
    {
        const Tag* p_;
    public:
        explicit TagIterator(const Tag &tag)
          : p_(&tag)
        {
        }

        explicit TagIterator(const unsigned char *tag)
          : p_(reinterpret_cast<const Tag *>(tag))
        {
        }

    private:
        friend class boost::iterator_core_access;

        void increment() { p_ = p_->end(); }

        bool equal(TagIterator const& other) const
        {
            return this->p_ == other.p_;
        }

        const Tag &dereference() const {
            return *p_; }
    };

    std::size_t blockLength() const
    {
        return common::extractLittleEndian<int32_t>(
                reinterpret_cast<const unsigned char *>(&block_size));
    }
    TagIterator tagsBegin() const
    {
        return TagIterator(getQual() + getLSeq());
    }
    TagIterator tagsEnd() const
    {
        return TagIterator(reinterpret_cast<const unsigned char *>(&refID) + blockLength());
    }
};

class BamParser
{
    unsigned headerBytesToSkip_;
    int referenceSequencesToSkip_;
public:
    BamParser() :
        headerBytesToSkip_(-1U),
        referenceSequencesToSkip_(-1) {}

    void reset()
    {
        headerBytesToSkip_ = -1U;
        referenceSequencesToSkip_ = -1;
    }

    template <typename CollectorT>
    bool parse(
        std::vector<char>::const_iterator &uncompressedIt,
        std::vector<char>::const_iterator uncompressedEnd,
        CollectorT &collector)
    {
        bool moreDataNeeded = true;
        if (headerBytesToSkip_)
        {
            skipHeader(uncompressedIt, uncompressedEnd);
        }

        if (!headerBytesToSkip_)
        {
            if (referenceSequencesToSkip_)
            {
                skipReferences(uncompressedIt, uncompressedEnd);
            }

            if (!referenceSequencesToSkip_)
            {
                while(uncompressedEnd != uncompressedIt)
                {
                    std::vector<char>::const_iterator last = uncompressedIt;
                    moreDataNeeded = parseBamRecord(uncompressedIt, uncompressedEnd, collector);
                    if(!moreDataNeeded || last == uncompressedIt)
                    {
                        // if the uncompressedIt did not move, means we can't parse from this point.
                        // if moreDataNeeded is false, the collector cannot accept anymore data
                        break;
                    }
                }
            }
        }

        return moreDataNeeded;
    }


private:
    bool skipHeader(
        std::vector<char>::const_iterator &it,
        std::vector<char>::const_iterator end);

    bool skipReferences(
        std::vector<char>::const_iterator &it,
        std::vector<char>::const_iterator end);

    template <typename ProcessorT>
    bool parseBamRecord(
        std::vector<char>::const_iterator &it,
        std::vector<char>::const_iterator end,
        ProcessorT &process)
    {
        static const unsigned BLOCK_SIZE_WIDTH = 4;
        const std::vector<char>::const_iterator blockIt = it;
        int block_size = 0;
        if (std::size_t(std::distance(it, end)) < BLOCK_SIZE_WIDTH)
        {
            return true;
        }

        it = common::extractLittleEndian(it, block_size);

        if (std::distance(it, end) < block_size)
        {
            it = blockIt;
            return true;
        }

        ISAAC_ASSERT_MSG(block_size >= int(sizeof(BamBlockHeader)), "bam record size is smaller than the minimum required block_size:" << block_size << " sizeof(BamBlockHeader):" << sizeof(BamBlockHeader));
        const BamBlockHeader &block = *reinterpret_cast<const BamBlockHeader *>(&*it - BLOCK_SIZE_WIDTH);

        const bool lastBlock = std::size_t(std::distance(it + block_size, end)) <= BLOCK_SIZE_WIDTH ||
            std::distance(it + block_size + BLOCK_SIZE_WIDTH, end) < (common::extractLittleEndian<unsigned>(it + block_size));
        const bool ret = process(block, lastBlock);
        it += block_size;

        return ret;
    }

};

/**
 * \brief extracts up to corresponding read length bcl sequence. if bam sequence is shorter than read length, the
 *        rest is padded with 0
 */
template <typename InsertIt, typename UnaryFun>
InsertIt extractBcl(
    const BamBlockHeader &bamBlock,
    InsertIt insertIt,
    UnaryFun translate,
    const flowcell::ReadMetadata &readMetadata)
{
    const unsigned char * const pSeq = bamBlock.getSeq();
    const unsigned char *pQual = bamBlock.getQual();
    const unsigned char * const pQualEnd = pQual + bamBlock.getLSeq();
    unsigned currentCycle = readMetadata.getFirstReadCycle();
    std::vector<unsigned>::const_iterator cycleIterator = readMetadata.getCycles().begin();

    unsigned seqOffset = 0;
    while(pQualEnd != pQual && readMetadata.getCycles().end() != cycleIterator)
    {
        if (*cycleIterator == currentCycle)
        {
            *insertIt++ = translate(bamToBcl(*pQual,
                                             (*(pSeq + seqOffset / 2) >> (4 * ((seqOffset + 1) % 2))) & 0x0F
                                             ));
            ++cycleIterator;
        }
        ++pQual;
        ++seqOffset;
        ++currentCycle;
    }
    insertIt = std::fill_n(insertIt, std::distance(cycleIterator, readMetadata.getCycles().end()), 0);
    return insertIt;
}

template <typename InsertIt>
InsertIt extractForwardBcl(
    const BamBlockHeader &bamBlock,
    InsertIt insertIt,
    const flowcell::ReadMetadata &readMetadata)
{
    return extractBcl(bamBlock, insertIt, [](unsigned char uc){return uc;}, readMetadata);
}

template <typename RandomAccessIt>
RandomAccessIt extractReverseBcl(
    const BamBlockHeader &bamBlock,
    RandomAccessIt randomAccessIt,
    const flowcell::ReadMetadata &readMetadata)
{
    extractBcl(
        bamBlock,
        boost::make_reverse_iterator(randomAccessIt + readMetadata.getLength()),
        &oligo::getReverseBcl,
        readMetadata).base();
    return randomAccessIt + readMetadata.getLength();
}

template <typename RandomAccessIt>
RandomAccessIt extractBcl(
    const BamBlockHeader &bamBlock,
    RandomAccessIt randomAccessIt,
    const flowcell::ReadMetadata &readMetadata)
{
    return bamBlock.isReverse() ?
        extractReverseBcl(bamBlock, randomAccessIt, readMetadata) :
        extractForwardBcl(bamBlock, randomAccessIt, readMetadata);
}


/**
 * \brief retrieve name of the read. Pad up to nameLengthMax with 0
 *
 * \return it + nameLengthMax
 */
template <typename ForwardAccessIt, typename InsertIt>
InsertIt extractReadName(
    ForwardAccessIt nameBegin, const std::size_t nameLength, const unsigned nameLengthMax, InsertIt insertIt)
{
    if (nameLengthMax > nameLength)
    {
        insertIt = std::copy(nameBegin, nameBegin + nameLength, insertIt);
        insertIt = std::fill_n(insertIt, nameLengthMax - nameLength, 0);
        //ISAAC_THREAD_CERR << "Extracting full:" <<  common::makeFastIoString(nameBegin, nameBegin + nameLength) << std::endl;
    }
    else
    {
        insertIt = std::copy(nameBegin + nameLength - nameLengthMax, nameBegin + nameLength, insertIt);
        //ISAAC_THREAD_CERR << "Extracting short:" <<  common::makeFastIoString(nameBegin + nameLength - nameLengthMax, nameBegin + nameLength) << std::endl;
    }
    return insertIt;
}

inline std::size_t formatReadName(
    const BamBlockHeader& block,
    char *szBuffer,
    const unsigned nameLengthMax,
    const flowcell::ReadMetadata &readMetadata)
{
    common::StaticVector<char, 1020> cigarBuffer;

    int pos = -1;
    if (readMetadata.getLength() != unsigned(block.getLSeq()))
    {
        // if input bam was --use-bases-masked, just make the CIGAR up.
        cigarBuffer.resize(100);
        snprintf(&cigarBuffer.front(), cigarBuffer.size(), "%dM", readMetadata.getLength());
    }
    else
    {
        for (BamBlockHeader::TagIterator it = block.tagsBegin(); block.tagsEnd() != it; ++it)
        {
            if (it->isTag("OC"))
            {
                std::copy(it->valueBegin(), it->valueEnd(), std::back_inserter(cigarBuffer));
            }
            else if (it->isTag("OP"))
            {
                // tag position is 1-based
                pos = it->value.i - 1;
            }
        }

        if (cigarBuffer.empty() && block.getCigarLength())
        {
            alignment::Cigar::toString(block.getCigar(), block.getCigar() + block.getCigarLength(), cigarBuffer);
            cigarBuffer.push_back(0);
        }

    }
    snprintf(
        szBuffer, nameLengthMax, "%c:%02d:%010d:%s",
        block.isUnmapped() ? 'u' : block.isReverse() ? 'r' : 'f',
        block.getRefId(), -1 == pos ? block.getPos() : pos,
        cigarBuffer.empty() ? "" : &cigarBuffer.front());

    return strlen(szBuffer);
}

/**
 * \brief retrieve name of the read. Pad up to nameLengthMax with 0
 *
 * \return it + nameLengthMax
 */
template <typename InsertIt>
InsertIt extractReadName(
    const BamBlockHeader &r1Block,
    const BamBlockHeader &r2Block,
    const unsigned nameLengthMax,
    const flowcell::ReadMetadata &r1Metadata,
    const flowcell::ReadMetadata &r2Metadata,
    InsertIt insertIt)
{
#ifdef ISAAC_DEV_STATS_ENABLED
    // ignore part of the name that follows the hash sign
    const char* hashPos = std::find(r1Block.nameBegin(), r1Block.nameEnd() - 1, '#');
    unsigned len = std::distance(r1Block.nameBegin(), hashPos) + 1;
    char szBuffer[std::max(nameLengthMax, len)];
    *std::copy(r1Block.nameBegin(), hashPos, szBuffer) = '#';

    if (len < nameLengthMax)
    {
        len += formatReadName(r1Block, szBuffer + len, nameLengthMax - len, r1Metadata);
        if (len < nameLengthMax)
        {
            szBuffer[len] = '-';
            ++len;
            len += formatReadName(r2Block, szBuffer + len, nameLengthMax - len, r2Metadata);
        }
    }

    return extractReadName(szBuffer, std::min(nameLengthMax, len), nameLengthMax, insertIt);
#else
    return extractReadName(r1Block.nameBegin(), r1Block.getReadNameLength(), nameLengthMax, insertIt);
#endif
}

template <typename InsertIt>
InsertIt extractReadName(
    const BamBlockHeader &bamBlock,
    const unsigned nameLengthMax,
    const flowcell::ReadMetadata &readMetadata,
    InsertIt insertIt)
{
#ifdef ISAAC_DEV_STATS_ENABLED
    ISAAC_ASSERT_MSG(!bamBlock.isPaired(), "TODO: support pairing of unpaired pairs");
    char szBuffer[nameLengthMax];
    const unsigned len = formatReadName(bamBlock, szBuffer, nameLengthMax, readMetadata);
    return extractReadName(szBuffer, std::min(nameLengthMax, len), nameLengthMax, insertIt);
#else
    return extractReadName(bamBlock.nameBegin(), bamBlock.getReadNameLength(), nameLengthMax, insertIt);
#endif
}


} // namespace bam
} // namespace isaac

#endif // #ifndef iSAAC_BAM_BAM_PARSER_HH
