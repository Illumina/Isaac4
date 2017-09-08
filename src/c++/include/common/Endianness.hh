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
 ** \file FastIo.hh
 **
 ** \brief Fast IO routines for integers and fixed width floating points.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_COMMON_ENDIANNESS_HH
#define iSAAC_COMMON_ENDIANNESS_HH


namespace isaac
{
namespace common
{

template <typename T, typename IterT> IterT extractLittleEndian(const IterT p, T &result)
{
#ifdef ISAAC_IS_BIG_ENDIAN
    // TODO:
    ISAAC_ASSERT_MSG(false, "TODO: implement little endian conversion");
#else
    result = *reinterpret_cast<const T*>(&*p);
#endif
    return p + sizeof(T);
}

template <typename T, typename IterT> T extractLittleEndian(const IterT p)
{
#ifdef ISAAC_IS_BIG_ENDIAN
    // TODO:
    ISAAC_ASSERT_MSG(false, "TODO: implement little endian conversion");
#endif
    return *reinterpret_cast<const T*>(&*p);
}

} // namespace common
} // namespace isaac

#endif // #ifndef iSAAC_COMMON_ENDIANNESS_HH
