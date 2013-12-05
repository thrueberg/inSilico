cmake_minimum_required (VERSION 2.6)
project (MixedPoisson)

set(CMAKE_BUILD_TYPE "Debug") # Release

# find the config directory
if(DEFINED ENV{INSILICOROOT})
  set(CONFDIR "$ENV{INSILICOROOT}/config")
else()
  set(CONFDIR "../../../config")
endif()

include( ${CONFDIR}/insilico.cmake )



# add the executable (name  sourceFile)
add_executable(mixedPoisson ../mixedPoisson.cpp)