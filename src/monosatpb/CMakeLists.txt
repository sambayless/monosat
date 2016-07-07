cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(minisatp)

#--------------------------------------------------------------------------------------------------
# Configurable options:

option(STATIC_BINARIES "Link binaries statically." ON)
option(USE_SORELEASE   "Use SORELEASE in shared library filename." ON)

#--------------------------------------------------------------------------------------------------
# Library version:

set(MINISATP_SOMAJOR   1)
set(MINISATP_SOMINOR   1)
set(MINISATP_SORELEASE 0)

# Compute VERSION and SOVERSION:
if(USE_SORELEASE)
  set(MINISATP_VERSION ${MINISATP_SOMAJOR}.${MINISATP_SOMINOR}.${MINISATP_SORELEASE})
else()
  set(MINISATP_VERSION ${MINISATP_SOMAJOR}.${MINISATP_SOMINOR})
endif()
set(MINISATP_SOVERSION ${MINISATP_SOMAJOR})

#--------------------------------------------------------------------------------------------------
# Compile flags:

add_definitions(-include Global.h -include Main.h -D_FILE_OFFSET_BITS=64)

#--------------------------------------------------------------------------------------------------
# Dependencies:

find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})

include_directories(${minisat_SOURCE_DIR})
include_directories(${minisatp_SOURCE_DIR})
include_directories(${minisatp_SOURCE_DIR}/ADTs)

#--------------------------------------------------------------------------------------------------
# Build Targets:

set(MINISATP_LIB_SOURCES
    ADTs/FEnv.cc
    ADTs/File.cc
    ADTs/Global.cc
    Debug.cc
    Hardware_adders.cc
    Hardware_clausify.cc
    Hardware_sorters.cc
    Main.cc
    PbParser.cc
    PbSolver.cc
    PbSolver_convert.cc
    PbSolver_convertAdd.cc
    PbSolver_convertBdd.cc
    PbSolver_convertSort.cc)

add_library(minisatp-lib-static STATIC ${MINISATP_LIB_SOURCES})
add_library(minisatp-lib-shared SHARED ${MINISATP_LIB_SOURCES})

target_link_libraries(minisatp-lib-shared minisat-lib-shared ${GMP_LIBRARY})
target_link_libraries(minisatp-lib-static minisat-lib-static ${GMP_LIBRARY})

set_target_properties(minisatp-lib-static PROPERTIES OUTPUT_NAME "minisatp")
set_target_properties(minisatp-lib-shared 
  PROPERTIES
    OUTPUT_NAME "minisatp" 
    VERSION ${MINISATP_VERSION}
    SOVERSION ${MINISATP_SOVERSION})

add_executable(minisatp Main.cc)

target_link_libraries(minisatp minisatp-lib-static)

if(STATIC_BINARIES)
  target_link_libraries(minisatp minisatp-lib-static)
else()
  target_link_libraries(minisatp minisatp-lib-shared)
endif()


#--------------------------------------------------------------------------------------------------
# Installation targets:
