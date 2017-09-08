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
 ** \file ContigLoader.cpp
 **
 ** Helper utility for loading multiple contigs of a fasta file.
 **
 ** \author Roman Petrovski
 **/
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "reference/ContigLoader.hh"

namespace isaac
{
namespace reference
{

void loadContig(
    const reference::SortedReferenceMetadata::Contig &contigMetadata,
    ContigList::UpdateRange &contig)
{
    ISAAC_ASSERT_MSG(contig.getLength() == contigMetadata.totalBases_, "Attempt to load wrong data into contig:" << contigMetadata << " " << contig)
    std::ifstream is(contigMetadata.filePath_.string().c_str());
    if (!is) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open reference file " + contigMetadata.filePath_.string()));
    }
    if (!is.seekg(contigMetadata.offset_))
    {
        using common::IoException;
        using boost::format;
        const format message = (boost::format("Failed to reach offset %d in reference file % s") % contigMetadata.offset_ % contigMetadata.filePath_);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }
//        ISAAC_THREAD_CERR << (boost::format("Contig seek %s (%3d:%8d): %s") % contigMetadata.name_ % contigMetadata.index_ % contigMetadata.totalBases_ % contigMetadata.filePath_).str() << std::endl;
    static const oligo::Translator<true> translator = {};
    char base = 0;
    std::size_t acgtBases = 0;
    auto updateIterator = contig.begin();
    while(is && (std::size_t(std::distance(contig.begin(), updateIterator)) < contigMetadata.totalBases_) && is.get(base))
    {
        if ('\r' != base && '\n' != base)
        {
            ISAAC_ASSERT_MSG(std::isalpha(base), "Invalid base read from " << contigMetadata << " : " << base);
            {
                ISAAC_ASSERT_MSG(updateIterator < contig.end(), "Trying to update contig past its range " << contig << " " << contigMetadata);
                *updateIterator = oligo::getBase(translator[base], true);
                acgtBases += oligo::REFERENCE_OLIGO_N != *updateIterator;
                ++updateIterator;
            }
        }
    }
    if (contigMetadata.totalBases_ != std::size_t(std::distance(contig.begin(), updateIterator)))
    {
        using common::IoException;
        using boost::format;
        const format message = (format("Failed to read %d bases from reference file % s: %d") % contigMetadata.totalBases_ % contigMetadata.filePath_.string() % contig.getLength());
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    if (contigMetadata.acgtBases_ != acgtBases)
    {
        using common::IoException;
        using boost::format;
        const format message = (format("Failed to read %d ACGT bases from reference file % s: %d") %
            contigMetadata.acgtBases_ % contigMetadata.filePath_.string() % acgtBases);
        BOOST_THROW_EXCEPTION(IoException(errno, message.str()));
    }

    //ISAAC_TRACE_STAT("Loaded contig " << contigMetadata << " ");
}

struct DummyFilter {bool operator() (const unsigned contigIdx) const {return true;}} dummyFilter;
/**
 * \brief loads all the fasta file contigs into memory on multiple threads
 */

} // namespace reference
} // namespace isaac
