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
 ** \file ReferenceMetadata.hh
 **
 ** Packaging of the metadata associated to a sorted reference.
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_FLOWCELL_REFERENCE_METADATA_HH
#define iSAAC_FLOWCELL_REFERENCE_METADATA_HH

#include <iostream>
#include <vector>

namespace isaac
{
namespace reference
{

/**
 ** \brief Read-only interface to the metadata associated to a reference.
 **
 ** The intended usage is for barcode management in ordered collections (the index
 ** in the collectionis is associated to each barcode metadata instance).
 ** Index 0 is reserved for mapping the barcode sequences that don't match the known ones
 **
 **/
class ReferenceMetadata
{
public:
    ReferenceMetadata(const std::string &name,
                      const boost::filesystem::path &path,
                      const unsigned index) :
        name_(name), path_(path), index_(index)
    {}

    const std::string &getName() const {return name_;}
    const boost::filesystem::path &getPath() const {return path_;}
    unsigned getIndex() const {return index_;}
    bool isXml() const
    {
        std::string ext = path_.extension().string();
        boost::to_upper(ext);
        return ext == ".XML";
    }

private:
    std::string name_;
    boost::filesystem::path path_;
    unsigned index_;
};

typedef std::vector<ReferenceMetadata> ReferenceMetadataList;

inline std::ostream &operator<<(std::ostream &os, const ReferenceMetadata &referenceMetadata)
{
    return os << "ReferenceMetadata(" << referenceMetadata.getName() << "," <<
        referenceMetadata.getPath() << "," << referenceMetadata.getIndex() << ")";
}

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_FLOWCELL_REFERENCE_METADATA_HH

