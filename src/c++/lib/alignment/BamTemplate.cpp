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
//  , debugClass_ = 0;
{
}

BamTemplate::BamTemplate(
    const FragmentMetadata &oneRead,
    const FragmentMetadata &anotherRead,
    const bool properPair/* = false*/,
    const unsigned alignmentScore/* = -1U*/)
    : alignmentScore_(alignmentScore)
    , pairInfo_(oneRead, anotherRead, properPair)
//  , debugClass_ = 0;
{
    if (!oneRead.getReadIndex())
    {
        fragmentMetadataList_.push_back(oneRead);
        fragmentMetadataList_.push_back(anotherRead);
    }
    else
    {
        fragmentMetadataList_.push_back(anotherRead);
        fragmentMetadataList_.push_back(oneRead);
    }
    ISAAC_ASSERT_MSG(oneRead.getReadIndex() != anotherRead.getReadIndex(), "Read index can't be same:" << oneRead << " " << anotherRead);
}

BamTemplate::BamTemplate(
    const FragmentMetadata &oneRead)
    : alignmentScore_(-1U)
//  , debugClass_ = 0;
{
    fragmentMetadataList_.push_back(oneRead);
}

BamTemplate::BamTemplate(const flowcell::ReadMetadataList &tileReads, const Cluster &cluster)
    : alignmentScore_(-1U)
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
        const FragmentMetadata &fragment = getFragmentMetadata(i);
        if (fragment.isAligned() && mapqThreshold <= fragment.mapQ)
        {
            ++goodFragments;
        }
    }

    return getFragmentCount() != goodFragments;
}

} // namespace alignment
} // namespace isaac
