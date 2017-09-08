################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2017 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## GNU GENERAL PUBLIC LICENSE Version 3
##
## You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
## along with this program. If not, see
## <https://github.com/illumina/licenses/>.
##
################################################################################
##
## file gcc.cmake
##
## CMake configuration file to configure global compiler settings
##
## author Roman Petrovski
##
################################################################################

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^x86_64$")
  if(NOT iSAAC_AVX2)
    set(iSAAC_VECTORIZATION "-msse4.2")
  else(NOT iSAAC_AVX2)
    set(iSAAC_VECTORIZATION "-mavx2")
  endif(NOT iSAAC_AVX2)
endif (CMAKE_SYSTEM_PROCESSOR MATCHES "^x86_64$")

set (iSAAC_CXX_OPTIMIZATION_FLAGS "${iSAAC_VECTORIZATION} -O2 -ftree-vectorize -finline-functions -fpredictive-commoning -fgcse-after-reload -funswitch-loops -ftree-slp-vectorize -fvect-cost-model -fipa-cp-clone -ftree-phiprop")
# -std=gnu++0x must be used to ensure cygwin builds succeed
set (CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} $ENV{CXXFLAGS} -std=gnu++0x -fpermissive -fopenmp -Wall -Wextra -Wunused -Wno-long-long -Wsign-compare -Wpointer-arith -DBOOST_SYSTEM_API_CONFIG_HPP -DBOOST_POSIX_API " CACHE STRING "g++ flags" FORCE)
# -03 causes loop unrolling that prevent autovectorization of some parts of BandedSmithWaterman

