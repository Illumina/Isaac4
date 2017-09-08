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
 ** \file DebugFragmentStorage.cpp
 **
 ** \author Roman Petrovski
 **/

#include <algorithm>
#include <atomic>

#include "alignment/matchSelector/debugStorage/MapqStatistics.hh"
#include "common/Debug.hh"
#include "common/FastIo.hh"

namespace isaac
{
namespace alignment
{
namespace matchSelector
{
namespace debugStorage
{

std::pair<BclClusters::const_iterator, BclClusters::const_iterator> getDebugName(const std::size_t readIndex, const FragmentMetadata &fragment)
{
    BclClusters::const_iterator debugNameStart = std::find(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd(), '#');
    debugNameStart = fragment.getCluster().nameEnd() != debugNameStart ? debugNameStart + 1 : fragment.getCluster().nameBegin();
    std::pair<BclClusters::const_iterator, BclClusters::const_iterator> ret = {debugNameStart, fragment.getCluster().nameEnd()};
    if (ret.second == ret.first)
    {
        return ret;
    }

    if (!readIndex)
    {
        ret.first = debugNameStart;
        ret.second = std::find(debugNameStart, fragment.getCluster().nameEnd(), '-');
    }
    else
    {
        ret.first = std::find(debugNameStart, fragment.getCluster().nameEnd(), '-') + 1;
        ret.second = fragment.getCluster().nameEnd();
    }

    return ret;
}

reference::ReferencePosition getAlignmentPositionFromName(const std::size_t readNumber, const FragmentMetadata &fragment)
{
    // numbers are 1-based
    const auto name = getDebugName(readNumber - 1, fragment);

    if (name.second == name.first)
    {
        return reference::ReferencePosition(reference::ReferencePosition::TooManyMatch);
    }

    const auto mateName = getDebugName(!(readNumber - 1), fragment);

    if (mateName.second == mateName.first)
    {
        return reference::ReferencePosition(reference::ReferencePosition::TooManyMatch);
    }

    // simulated data is assumed to be paired for the moment. With real data, assume that input singletons are aligned correctly since
    // they are most likely incorrectly aligned in the input anyway.
    if ('u' == *name.first || 'u' == *mateName.first)
    {
//        ISAAC_ASSERT_MSG(false, common::makeFastIoString(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()) << " " << fragment);
//        return reference::ReferencePosition(reference::ReferencePosition::NoMatch);
        // unaligned are present in bams from real aligners such as bwa. Pretend that whatever they don't align we align correctly
        return reference::ReferencePosition(reference::ReferencePosition::TooManyMatch);
    }
    return reference::ReferencePosition(
        std::atol(&*name.first + 2),
        std::atol(&*std::find(name.first + 2, name.second, ':') + 1),
		'r' == *name.first);
}

std::pair<const char *, const char *> getCigarFromName(const std::size_t readNumber, const FragmentMetadata &fragment)
{
    // numbers are 1-based
    const auto name = getDebugName(readNumber - 1, fragment);

    if (name.second == name.first)
    {
        return std::make_pair<const char *, const char *>(0,0);
    }

    if ('u' == *name.first)
    {
        return std::make_pair<const char *, const char *>(0,0);
    }

    auto ret = std::make_pair<const char *, const char *>(&*name.first, &*name.second);

    ret.first = std::find(ret.first, ret.second, ':');
    if (ret.second != ret.first)
    {
        ret.first = std::find(ret.first + 1, ret.second, ':');
        if (ret.second != ret.first)
        {
            ret.first = std::find(ret.first + 1, ret.second, ':');
            ret.first += (ret.second != ret.first);
        }
    }
    return ret;
}


bool compareCigars(const std::size_t readNumber, const FragmentMetadata &fragment)
{
    common::StaticVector<char, 1020> cigarBuffer;
    alignment::Cigar::toString(fragment.cigarBegin(), fragment.cigarEnd(), cigarBuffer);

    std::pair<const char *, const char *> nameCigar = getCigarFromName(readNumber, fragment);

//    ISAAC_THREAD_CERR << " name " << common::makeFastIoString(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()) << std::endl;
//    ISAAC_THREAD_CERR << " nameCigar " << common::makeFastIoString(nameCigar.first, nameCigar.second) << std::endl;
//    ISAAC_THREAD_CERR << " cigarBuffer " << common::makeFastIoString(cigarBuffer.begin(), cigarBuffer.end()) << std::endl;
    if (cigarBuffer.size() == std::size_t(std::distance(nameCigar.first, nameCigar.second)))
    {
        if (nameCigar.second == std::mismatch(nameCigar.first, nameCigar.second, cigarBuffer.begin()).first)
        {
            return true;
        }
    }
    return false;
}

bool alignsCorrectly(const std::size_t readNumber, const FragmentMetadata &fragment)
{
    if (alignment::containsHomopolymer(fragment.getStrandSequence().begin(), fragment.getStrandSequence().end()))
    {
        return true;
    }
    const reference::ReferencePosition oriPos = getAlignmentPositionFromName(readNumber, fragment);
    if (oriPos.isTooManyMatch())
    {
        return true;
    }
//    ISAAC_THREAD_CERR << "oriPos:" << oriPos << " name " << common::makeFastIoString(fragment.getCluster().nameBegin(), fragment.getCluster().nameEnd()) << std::endl;
    return fragment.isReverse() == oriPos.reverse() &&
        fragment.getContigId() == oriPos.getContigId() &&
        uint64_t(fragment.getPosition()) == oriPos.getPosition() &&
        compareCigars(readNumber, fragment);
}

} // namespace debugStroage
} // namespace matchSelector
} // namespace alignment
} // namespace isaac
