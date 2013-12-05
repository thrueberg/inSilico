# ------------------------------------------------------------------------------
# find root directory via variable or set locally
if(DEFINED ENV{INSILICOROOT})
  set(ROOTDIR "$ENV{INSILICOROOT}")
else()
  set(ROOTDIR "..")
endif()

message(STATUS "Using Insilico root directory: " "${ROOTDIR}")
include_directories(${ROOTDIR})

# ------------------------------------------------------------------------------
# find Eigen3
set(CMAKE_MODULE_PATH ${ROOTDIR}/config)  # pass to other directory
find_package(Eigen3 REQUIRED)
message( STATUS "Eigen3 version: " ${EIGEN3_VERSION} )
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# ------------------------------------------------------------------------------
# find Boost
find_package(Boost REQUIRED)
include_directories(SYSTEM ${BOOST_INCLUDE_DIR})


# ------------------------------------------------------------------------------
# check compilation environment
message( STATUS "CMake compiler ID: " ${CMAKE_CXX_COMPILER_ID})

if (     "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # ----------------------------------------------------------------------------
  # using GCC
  message( STATUS "Using GNU compiler" )

  # standard flags
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
  # debug flags
  set ( ADDITIONAL_DEBUG_FLAGS
    "-g3 -Wextra -ansi -pedantic -Wno-unused-parameter -Wconversion -std=c++0x")
  set ( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG ${ADDITIONAL_DEBUG_FLAGS}" )
  # release flags
  set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp" )
  
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # ----------------------------------------------------------------------------
  # using Intel C++
  
  # Currently invoked by: 
  # cmake -D CMAKE_C_COMPILER=/opt/intel/bin/icc -D CMAKE_CXX_COMPILER=/opt/intel/bin/icpc .
  # or:
  # CC=/opt/intel/bin/icc CXX=/opt/intel/bin/icpc cmake .
  # make sure to 'purge' first

  message( STATUS "Using GNU compiler" )

  # standard flags
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
  # debug flags
  set ( ADDITIONAL_DEBUG_FLAGS "-ansi -wd383,869,981,1418,1419,2196 -mkl=sequential" )
  set ( CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG ${ADDITIONAL_DEBUG_FLAGS}" )
  # release flags
  set ( ADDITIONAL_RELEASE_FLAGS "-openmp -wd858 -mkl=parallel" )
  set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_RELEASE_FLAGS}" )

elseif ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # ----------------------------------------------------------------------------
  # using Clang
  message( FATAL_ERROR "Clang++ is not yet supported" )

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # ----------------------------------------------------------------------------
  # using Visual Studio C++
  message( FATAL_ERROR "MSVC is never supported" )

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MinGW")
  # ----------------------------------------------------------------------------
  # using MinGW C++
  message( FATAL_ERROR "MinGW is not yet supported" )

  # set flags
endif()



# ------------------------------------------------------------------------------
# Additional target removes the files generate by CMake
add_custom_target (
    purge
    COMMAND rm -rf CMakeFiles && rm CMakeCache.txt cmake_install.cmake Makefile
 )