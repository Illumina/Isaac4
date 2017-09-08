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
 ** \file SortedReferenceMetadata.cpp
 **
 ** Information about the pre-processed reference data files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "reference/SortedReferenceMetadata.hh"

namespace isaac
{
namespace reference
{

const unsigned SortedReferenceMetadata::OLDEST_SUPPORTED_REFERENCE_FORMAT_VERSION;
const unsigned SortedReferenceMetadata::CURRENT_REFERENCE_FORMAT_VERSION;

void SortedReferenceMetadata::makeAbsolutePaths(const boost::filesystem::path &basePath)
{
    BOOST_FOREACH(AllMaskFiles::value_type &maskFiles, maskFiles_)
    {
        BOOST_FOREACH(MaskFile &maskFile, maskFiles.second)
        {
            maskFile.path = boost::filesystem::absolute(maskFile.path, basePath);
        }
    }

    BOOST_FOREACH(AnnotationFile &annotationFile, annotationFiles_)
    {
        annotationFile.path_ = boost::filesystem::absolute(annotationFile.path_, basePath);
    }

    BOOST_FOREACH(Contig &contig, contigs_)
    {
        contig.filePath_ = boost::filesystem::absolute(contig.filePath_, basePath);
    }
}

void SortedReferenceMetadata::putContig(
    const Contig &contig)
{
    contigs_.push_back(contig);
}

void SortedReferenceMetadata::putContig(
    const uint64_t genomicOffset,
    const std::string& name,
    const boost::filesystem::path &sequencePath,
    const uint64_t byteOffset,
    const uint64_t byteSize,
    const uint64_t totalBases,
    const uint64_t acgtBases,
    const unsigned index,
    const std::string &bamSqAs,
    const std::string &bamSqUr,
    const std::string &bamM5
    )
{
    contigs_.push_back(Contig(
        index,
        name,
        false,
        sequencePath,
        byteOffset,
        byteSize,
        genomicOffset,
        totalBases, acgtBases,
        bamSqAs,
        bamSqUr,
        bamM5));
}
void SortedReferenceMetadata::addMaskFile(
    const unsigned seedLength,
    const unsigned int maskWidth,
    const unsigned mask, const boost::filesystem::path &filePath,
    const std::size_t kmers)
{
    maskFiles_[seedLength].push_back(MaskFile(filePath, maskWidth, mask, kmers));
}

uint64_t SortedReferenceMetadata::getTotalKmers(const unsigned seedLength) const
{
    unsigned maskWidth = -1U;
    uint64_t ret = 0;
    BOOST_FOREACH(const MaskFile &mask, maskFiles_.at(seedLength))
    {
        if (-1U == maskWidth)
        {
            maskWidth = mask.maskWidth;
        }
        ISAAC_ASSERT_MSG(maskWidth == mask.maskWidth, "Mixed mask widths are not supported");
        ret += mask.kmers;
    }
    return ret;
}

void SortedReferenceMetadata::merge(SortedReferenceMetadata &that)
{
    ISAAC_ASSERT_MSG(formatVersion_ == that.formatVersion_, "Incompatible formats: " << formatVersion_ << " and " << that.formatVersion_ << " cannot be merged");
    if(!contigs_.empty() && !that.contigs_.empty())
    {
        if (contigs_ != that.contigs_)
        {
            BOOST_THROW_EXCEPTION(common::FeatureNotAvailable("Cannot merge references with different contig lists"));
        }
    }
    else
    {
        contigs_.insert(contigs_.end(), that.contigs_.begin(), that.contigs_.end());
    }

    BOOST_FOREACH(AllMaskFiles::value_type &seedMaskFiles, that.maskFiles_)
    {
        MaskFiles &maskFiles = getMaskFileList(seedMaskFiles.first);
        maskFiles.insert(maskFiles.end(), seedMaskFiles.second.begin(), seedMaskFiles.second.end());
    }

    annotationFiles_.insert(annotationFiles_.end(), that.annotationFiles_.begin(), that.annotationFiles_.end());
//    if (annotationFiles_.empty())
//    {
//        annotationFiles_ = that.annotationFiles_;
//    }
//    else
//    {
//        ISAAC_ASSERT_MSG(that.annotationFiles_.empty(), "Only one of the references is expected to have annotation files");
//    }
}

bool SortedReferenceMetadata::singleFileReference() const {
    ISAAC_ASSERT_MSG(0 != getContigsCount(), "Contigs list must not be empty.");
    Contigs::const_iterator different = std::find_if(contigs_.begin() + 1, contigs_.end(),
                                                     boost::bind(&Contig::filePath_, _1) != contigs_.front().filePath_);
    return contigs_.end() == different;
}

} // namespace reference
} // namespace isaac

