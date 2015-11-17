#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

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
# Thirdparty libraries
#------------------------------------------------------------------------------#

set( TPL_INSTALL_PREFIX /path/to/third/party/install 
                        CACHE PATH
                        "path to thirdparty install" )
if (NOT TPL_INSTALL_PREFIX STREQUAL "")
  set(EXODUS_ROOT  ${TPL_INSTALL_PREFIX})
endif()


find_package(Boost 1.47 REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )


#------------------------------------------------------------------------------#
# Enable IO with exodus
#------------------------------------------------------------------------------#

option(ENABLE_IO "Enable I/O with third party libraries." OFF)
if(ENABLE_IO)

  set( IO_LIBRARIES )

  find_library ( EXODUS_LIBRARY 
                 NAMES exodus 
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  find_path    ( EXODUS_INCLUDE_DIR 
                 NAMES exodusII.h
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES include
                 NO_DEFAULT_PATH )

  find_library ( NETCDF_LIBRARY 
                 NAMES netcdf 
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  find_library ( HDF5_LIBRARY 
                 NAMES hdf5 
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  find_library ( HDF5_HL_LIBRARY 
                 NAMES hdf5_hl 
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  find_library ( SZIP_LIBRARY 
                 NAMES szip 
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  find_library ( Z_LIBRARY 
                 NAMES z
                 PATHS ${EXODUS_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  if (EXODUS_LIBRARY AND EXODUS_INCLUDE_DIR) 
     message(STATUS "Found Exodus: ${EXODUS_ROOT}")
     set( EXODUS_FOUND TRUE )
     list( APPEND IO_LIBRARIES ${EXODUS_LIBRARY}
                               ${NETCDF_LIBRARY}
                               ${HDF5_HL_LIBRARY}
                               ${HDF5_LIBRARY}
                               ${SZIP_LIBRARY}
                               ${Z_LIBRARY}
                               -ldl )
     include_directories( ${EXODUS_INCLUDE_DIR} )
     add_definitions( -DHAVE_EXODUS )
  endif()

  if ( NOT IO_LIBRARIES )
     MESSAGE( FATAL_ERROR "Need to specify EXODUS" )
  endif()

endif(ENABLE_IO)


#~---------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
