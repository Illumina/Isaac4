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
 ** \file BarcodePathMap.hh
 **
 ** Helper class for mapping barcodes to output files.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BARCODE_BAM_HH
#define iSAAC_BUILD_BARCODE_BAM_HH

#include <iterator>
#include <boost/filesystem.hpp>

#include "common/FileSystem.hh"
#include "flowcell/BarcodeMetadata.hh"

namespace isaac
{
namespace demultiplexing
{

class BarcodePathMap
{
public:
    typedef std::vector<unsigned> BarcodeSampleIndexMap;
    typedef std::vector<unsigned> BarcodeProjectIndexMap;

    BarcodePathMap() : projectIndexMax_(-1U){}
    /**
     * \param projectIds one entry per barcode index mapping it to the corresponding project id
     * \param sampleIds one entry per barcode index mapping it to the corresponding sample
     * \param samplePaths one entry per sample id
     */
    BarcodePathMap(
        const BarcodeProjectIndexMap &projectIds,
        const BarcodeSampleIndexMap &sampleIds,
        const std::vector<boost::filesystem::path> &samplePaths):
            barcodeProjectIndex_(projectIds),
            projectIndexMax_(std::distance(barcodeProjectIndex_.begin(), std::max_element(barcodeProjectIndex_.begin(), barcodeProjectIndex_.end()))),
            barcodeSampleIndex_(sampleIds), samplePaths(samplePaths){}
    /// Each position in the vector contains unique index of the project-sample
    const BarcodeSampleIndexMap &getSampleIndexMap() const {return barcodeSampleIndex_;}
    const std::vector<boost::filesystem::path> &getPaths() const {return samplePaths;}
    unsigned getTotalBarcodes() const {return barcodeSampleIndex_.size();}
    unsigned getTotalSamples() const {return samplePaths.size();}
    unsigned getProjectIndex(const unsigned barcodeIndex) const {return barcodeProjectIndex_.at(barcodeIndex);}
    unsigned getMaxProjectIndex() const {return projectIndexMax_;}
    unsigned getSampleIndex(const unsigned barcodeIndex) const {return barcodeSampleIndex_.at(barcodeIndex);}
    const boost::filesystem::path &getSampleFilePath(std::size_t sampleId) const
    {
        return samplePaths.at(sampleId);
    }
    const boost::filesystem::path &getFilePath(const flowcell::BarcodeMetadata &barcode) const
    {
        return getSampleFilePath(getSampleIndex(barcode.getIndex()));
    }


private:
    template<class Archive> friend void serialize(Archive & ar, BarcodePathMap &, const unsigned int file_version);
    BarcodeProjectIndexMap barcodeProjectIndex_;
    unsigned projectIndexMax_;
    BarcodeSampleIndexMap barcodeSampleIndex_;
    std::vector<boost::filesystem::path> samplePaths;
};


inline boost::filesystem::path getBarcodeFilePath(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadata &barcode,
    const char *fileName)
{
    return outputDirectory / barcode.getProject() / barcode.getSampleName() / fileName;
}

/**
 * \brief Produces mapping so that all barcodes having the same sample name go into the same output file.
 *
 * \return Returns the pair of a vector of that maps a barcode index to a unique output file index in the
 *         second vector so that two barcodes that are supposed to go into the same file will end up
 *         having the same mapping
 */
inline BarcodePathMap mapBarcodesToFiles(
    const boost::filesystem::path &outputDirectory,
    const flowcell::BarcodeMetadataList &barcodeMetadataList,
    const char *fileName)
{
    // Map barcodes to projects
    std::vector<std::string> projects;
    std::transform(barcodeMetadataList.begin(), barcodeMetadataList.end(), std::back_inserter(projects),
                   boost::bind(&flowcell::BarcodeMetadata::getProject, _1));
    std::sort(projects.begin(), projects.end());
    projects.erase(std::unique(projects.begin(), projects.end()), projects.end());

    BarcodePathMap::BarcodeProjectIndexMap barcodeProject(barcodeMetadataList.size());
    for(const flowcell::BarcodeMetadata &barcode : barcodeMetadataList)
    {
        barcodeProject.at(barcode.getIndex()) =
            std::distance(projects.begin(), std::lower_bound(projects.begin(), projects.end(), barcode.getProject()));
    }

    // Map barcodes to sample paths
    std::vector<boost::filesystem::path> samples;
    std::transform(barcodeMetadataList.begin(), barcodeMetadataList.end(), std::back_inserter(samples),
                   boost::bind(&getBarcodeFilePath, outputDirectory, _1, fileName));
    std::sort(samples.begin(), samples.end());
    samples.erase(std::unique(samples.begin(), samples.end()), samples.end());

    BarcodePathMap::BarcodeProjectIndexMap barcodeSample(barcodeMetadataList.size());
    for (const flowcell::BarcodeMetadata &barcode : barcodeMetadataList)
    {
        barcodeSample.at(barcode.getIndex()) =
            std::distance(samples.begin(), std::lower_bound(samples.begin(), samples.end(),
                    getBarcodeFilePath(outputDirectory, barcode, fileName)));
    }

    return BarcodePathMap(barcodeProject, barcodeSample, samples);
}

inline void createDirectories(
        const BarcodePathMap &barcodeBamMapping,
        const flowcell::BarcodeMetadataList& barcodeMetadataList)
{
    std::vector<boost::filesystem::path> directories;
    for (const flowcell::BarcodeMetadata &barcode : barcodeMetadataList)
    {
        const boost::filesystem::path &bamPath = barcodeBamMapping.getFilePath(barcode);
        directories.push_back(bamPath.parent_path());
        directories.push_back(directories.back().parent_path());
    }

    common::createDirectories(directories);
}


} // namespace demultiplexing
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BARCODE_BAM_HH
