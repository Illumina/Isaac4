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
 ** \file BamTemplate.cpp
 **
 ** \brief See BamTemplate.hh
 ** 
 ** \author Come Raczy
 **/

#include <boost/foreach.hpp>

#include "alignment/BamTemplate.hh"

namespace isaac
{
namespace alignment
{

BamTemplate::BamTemplate()
    : alignmentScore_(-1U)
    , properPair_(false)
//  , debugClass_ = 0;
{
}

BamTemplate::BamTemplate(
    const FragmentMetadata &read1,
    const FragmentMetadata &read2,
    const bool properPair/* = false*/)
    : alignmentScore_(-1U)
    , properPair_(properPair)
//  , debugClass_ = 0;
{
    ISAAC_ASSERT_MSG(!read1.getReadIndex(), "Incorrect index for read1:" << read1);
    ISAAC_ASSERT_MSG(1 == read2.getReadIndex(), "Incorrect index for read2:" << read2);
    fragmentMetadataList_.push_back(read1);
    fragmentMetadataList_.push_back(read2);
}

BamTemplate::BamTemplate(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster)
    : alignmentScore_(-1U)
    , properPair_(false)
//  , debugClass_ = 0;
{
    BOOST_FOREACH(const flowcell::ReadMetadata &read, tileReads)
    {
        fragmentMetadataList_.push_back(FragmentMetadata(&cluster, read.getIndex()));
    }
}

void BamTemplate::reset(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster)
{
    *this = BamTemplate(tileReads, cluster);
}

bool BamTemplate::filterLowQualityFragments(const int mapqThreshold)
{
    if (-1 == mapqThreshold)
    {
        return false;
    }
    unsigned goodFragments = 0;
    for (unsigned i = 0; getFragmentCount() > i; ++i)
    {
        FragmentMetadata &fragment = getFragmentMetadata(i);
        if (fragment.isAligned() && mapqThreshold <= fragment.mapQ)
        {
            ++goodFragments;
        }
    }

    return getFragmentCount() != goodFragments;
}

} // namespace alignment
} // namespace isaac
