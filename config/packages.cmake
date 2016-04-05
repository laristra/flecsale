#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Thirdparty libraries
#------------------------------------------------------------------------------#

set( TPL_INSTALL_PREFIX /path/to/third/party/install 
                        CACHE PATH
                        "path to thirdparty install" )
if (TPL_INSTALL_PREFIX)
  set(EXODUS_ROOT ${TPL_INSTALL_PREFIX})
  set(TECIO_ROOT  ${TPL_INSTALL_PREFIX})
  set(ShaPo_DIR   ${TPL_INSTALL_PREFIX})
  set(VTK_DIR     ${TPL_INSTALL_PREFIX})
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

  find_library ( TECIO_LIBRARY 
                 NAMES tecio 
                 PATHS ${TECIO_ROOT} 
                 PATH_SUFFIXES lib
                 NO_DEFAULT_PATH )

  find_path    ( TECIO_INCLUDE_DIR 
                 NAMES TECIO.h
                 PATHS ${TECIO_ROOT} 
                 PATH_SUFFIXES include
                 NO_DEFAULT_PATH )

  if (TECIO_LIBRARY AND TECIO_INCLUDE_DIR) 
     message(STATUS "Found TECIO: ${TECIO_ROOT}")
     set( TECIO_FOUND TRUE )
     list( APPEND IO_LIBRARIES ${TECIO_LIBRARY} )
     include_directories( ${TECIO_INCLUDE_DIR} )
     add_definitions( -DHAVE_TECIO )
  endif()

  if ( NOT IO_LIBRARIES )
     MESSAGE( FATAL_ERROR "Need to specify EXODUS" )
  endif()

endif(ENABLE_IO)



#------------------------------------------------------------------------------#
# Enable shapo
#------------------------------------------------------------------------------#

option(ENABLE_VORO "Enable vornoi with shapo." OFF)
if(ENABLE_VORO)

  find_package(ShaPo HINTS ${ShaPo_DIR} )
  
  if (ShaPo_FOUND)
    message(STATUS "Found ShaPo: ${ShaPo_DIR}")
    include_directories( ${SHAPO_INCLUDE_DIRS} )  
    include( ${ShaPo_CONFIG} )
    include( ${SHAPO_USE_FILE} )
    add_definitions( -DHAVE_SHAPO )
  else ()
    message( FATAL_ERROR "Need to specify SHAPO" )
  endif()

  if (ShaPo_FOUND)
    find_package(VTK REQUIRED HINTS ${VTK_DIR} )
    if (VTK_FOUND)
      message(STATUS "Found VTK: ${VTK_DIR}")
      include( ${VTK_CONFIG} )
      include( ${VTK_USE_FILE} )
      add_definitions( -DHAVE_VTK )    
  else()
    message(STATUS "Error: VTK could not be found.")
  endif()

endif()

endif(ENABLE_VORO)


#~---------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
