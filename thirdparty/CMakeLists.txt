############################################################
cmake_minimum_required(VERSION 3.1.0)
############################################################

project(ThirdParty LANGUAGES CXX C Fortran)

# add cmake user module directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules/")

# other important directories
set( THIRDPARTY_PATCH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/patches/" )
set( THIRDPARTY_FILES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/files/" )

# installation path for third party libraries
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install"
    CACHE PATH "default install path" FORCE)
endif()

include(macros)


############################################################
# GENERAL OPTIONS
############################################################

include(CMakeDependentOption)

# third party libraries to build
option(ENABLE_FLECSI "Enable flecsi and dependencies" ON)


option(ENABLE_EXODUS "Enable Exodus and dependencies" OFF)
cmake_dependent_option(ENABLE_HDF5 "Enable HDF5" ON
                       "ENABLE_EXODUS" OFF)
cmake_dependent_option(ENABLE_NETCDF "Enable HDF5" ON
                       "ENABLE_EXODUS" OFF)
cmake_dependent_option(ENABLE_SZIP "Enable szip" ON
                       "ENABLE_HDF5" OFF)
cmake_dependent_option(ENABLE_ZLIB "Enable zlib" ON
                       "ENABLE_HDF5" OFF)

option(ENABLE_METIS "Enable METIS and dependencies" OFF)

option(ENABLE_SCOTCH "Enable SCOTCH and dependencies" OFF)

# cmake functionality
option(TPL_DOWNLOAD_ONLY "Download libraries to TPL_DOWNLOAD_PATH and exit" OFF)
set( TPL_DOWNLOAD_PATH "" CACHE FILEPATH "\
If TPL_DOWNLOAD_ONLY=ON, the location where files will be downloaded to. \
Otherwise, use thirdparty libraries that were previously downloaded here \ 
for the installation process.\
" )


############################################################
# FOR DOWNLOAD OPLY 
############################################################

if ( TPL_DOWNLOAD_ONLY )

  message(STATUS "Downloading thirdparty libraries")

  if ( NOT TPL_DOWNLOAD_PATH STREQUAL "" )
    set( DOWNLOAD_PATH_ ${TPL_DOWNLOAD_PATH} )
  else()
    message(WARNING "TPL_DOWNLOAD_PATH not set, using default." )
    set( DOWNLOAD_PATH_ ${THIRDPARTY_FILES_DIR} )
  endif ()
  message(STATUS "Downloading to: ${DOWNLOAD_PATH_}" )

  download_libraries( ${DOWNLOAD_PATH_} )

  return() # exit now

else()

  message(STATUS "Building thirdparty libraries")

endif()




############################################################
# BUILD LIBRARIES
############################################################

include(ExternalProject) 

include(versions)

include(zlib)
include(szip)
include(hdf5)
include(netcdf)
include(exodus)
include(metis)
include(scotch)

# 
# 
# if(ENABLE_FLECSI)
#   include(cmake/flecsi.cmake)
# endif()
