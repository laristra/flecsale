#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#

cinch_minimum_required(1.0)

#------------------------------------------------------------------------------#
# Set the project name
#------------------------------------------------------------------------------#

project(FleCSALE)

#------------------------------------------------------------------------------#
# Begin project setup
#------------------------------------------------------------------------------#

#if(COMMAND cmake_policy)
#  cmake_policy(SET CMP0005 NEW)  # generate escape sequences for defines
#  cmake_policy(SET CMP0012 NEW)  # recognize number & boolean literals
#endif(COMMAND cmake_policy)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# set some global variables
set( ALE_LIBRARIES )
set( ALE_DATA_DIR "${CMAKE_SOURCE_DIR}/data" )  
set( ALE_TOOL_DIR "${CMAKE_SOURCE_DIR}/tools" )  

#------------------------------------------------------------------------------#
# If a C++14 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#

include(cxx14)

check_for_cxx14_compiler(CXX14_COMPILER)

if(CXX14_COMPILER)
    enable_cxx14()
else()
    message(FATAL_ERROR "C++14 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Some precision setup
#------------------------------------------------------------------------------#

# double or single precision
OPTION (DOUBLE_PRECISION "Use double precision reals"  ON)

if( DOUBLE_PRECISION ) 
  message(STATUS "Note: Double precision build activated.")
  add_definitions( -DDOUBLE_PRECISION )
  SET (TEST_TOLERANCE 1.0e-14 CACHE STRING "The testing tolerance" )
else()
  message(STATUS "Note: Single precision build activated.")
  SET (TEST_TOLERANCE 1.0e-6 CACHE STRING "The testing tolerance" )
endif()

add_definitions( -DTEST_TOLERANCE=${TEST_TOLERANCE} )


# size of integer ids to use
option( USE_64BIT_IDS "Type of integer to use for ids" ON )

if( USE_64BIT_IDS ) 
  message(STATUS "Note: using 64 bit integer ids.")
  add_definitions( -DUSE_64BIT_IDS )
else()
  message(STATUS "Note: using 32 bit integer ids.")
endif()

#------------------------------------------------------------------------------#
# Enable Boost.Preprocessor
#------------------------------------------------------------------------------#

# This changes the Cinch default
set(ENABLE_BOOST_PREPROCESSOR ON CACHE BOOL
    "Enable Boost.Preprocessor")

#------------------------------------------------------------------------------#
# Add subprojects
#------------------------------------------------------------------------------#

set(FLECSI_RUNTIME_MODEL "serial" CACHE STRING
  "Select the runtime model")
set_property(CACHE FLECSI_RUNTIME_MODEL
  PROPERTY STRINGS serial mpilegion legion mpi)

find_package(FLECSI QUIET)
if(NOT FLECSI_FOUND)
  cinch_add_subproject("flecsi")
  list( APPEND ALE_LIBRARIES flecsi )
  include_directories( ${CMAKE_BINARY_DIR} )
else()
  include_directories(${FleCSI_INCLUDE_DIRS})
  list( APPEND ALE_LIBRARIES ${FleCSI_LIBRARIES})
endif()

#------------------------------------------------------------------------------#
# Add library targets
#------------------------------------------------------------------------------#

cinch_add_library_target(flecsale libsrc)
list( APPEND ALE_LIBRARIES flecsale )

#------------------------------------------------------------------------------#
# Set application directory
#------------------------------------------------------------------------------#

cinch_add_application_directory(apps)

#------------------------------------------------------------------------------#
# Set header suffix regular expression
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

