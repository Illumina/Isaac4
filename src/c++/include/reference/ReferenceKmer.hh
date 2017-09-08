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
 ** \file ReferenceKmer.hh
 **
 ** Representation of a k-mer at a given position in a reference genome.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
#define iSAAC_REFERENCE_REFERENCE_KMER_HH

#include <utility>

#include <boost/mpl/back.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/modulus.hpp>

#include "oligo/Kmer.hh"
#include "oligo/Permutate.hh"
#include "reference/ReferencePosition.hh"

namespace isaac
{
namespace reference
{

typedef boost::mpl::copy_if<
    oligo::SUPPORTED_KMERS
    , boost::mpl::equal_to<boost::mpl::modulus<boost::mpl::_1, boost::mpl::int_<4> >, boost::mpl::int_<0> >
    , boost::mpl::back_inserter< boost::mpl::vector<> >
    >::type SUPPORTED_KMERS;

static const unsigned FIRST_SUPPORTED_KMER = boost::mpl::front<SUPPORTED_KMERS>::type::value;
static const unsigned LAST_SUPPORTED_KMER = boost::mpl::back<SUPPORTED_KMERS>::type::value;

} // namespace reference
} // namespace isaac

#endif // #ifndef iSAAC_REFERENCE_REFERENCE_KMER_HH
