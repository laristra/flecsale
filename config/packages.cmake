#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#
#------------------------------------------------------------------------------#
# Thirdparty libraries
#------------------------------------------------------------------------------#

# find_package(Boost 1.47 REQUIRED)
# include_directories( ${Boost_INCLUDE_DIRS} )

set( IO_LIBRARIES )
set( VORO_LIBRARIES )

#------------------------------------------------------------------------------#
# Find libraries
#------------------------------------------------------------------------------#


find_package(EXODUSII)

find_package(TECIO QUIET)

find_package(ShaPo QUIET)
if (ShaPo_FOUND)
  message(STATUS "Found ShaPo: ${ShaPo_DIR}")
endif()

find_package(VTK QUIET)
if (VTK_FOUND)
  message(STATUS "Found VTK: ${VTK_DIR}")
endif()


# find python for running regression tests
FIND_PACKAGE (PythonInterp)
IF (PYTHONINTERP_FOUND)
   MESSAGE (STATUS "Found PythonInterp: ${PYTHON_EXECUTABLE}")
ELSE (PYTHONINTERP_FOUND)
   MESSAGE (FATAL "Did not find python. Python is needed to run regression tests.")
ENDIF ()

# find python for embedding
FIND_PACKAGE (PythonLibs)
IF (PYTHONLIBS_FOUND)
   MESSAGE (STATUS "Found PythonLibs: ${PYTHON_INCLUDE_DIRS}")
ELSE (PYTHONINTERP_FOUND)
   MESSAGE (FATAL "Did not find python. Python is needed to run regression tests.")
ENDIF ()


#------------------------------------------------------------------------------#
# Enable shapo
#------------------------------------------------------------------------------#

if ( ShaPo_FOUND AND VTK_FOUND )
  option(ENABLE_VORO "Enable vornoi with shapo." ON)
else()
  option(ENABLE_VORO "Enable vornoi with shapo." OFF)
endif()

if(ENABLE_VORO)
  
  if (VTK_FOUND)
    set( VTK_INCLUDED TRUE )
    include( ${VTK_CONFIG} )
    include( ${VTK_USE_FILE} )
    list(APPEND VORO_LIBRARIES ${VTK_LIBRARIES} )
    add_definitions( -DHAVE_VTK )    
  else()
    message( FATAL_ERROR "You need libVTK either from TPL or system to enable voro")
  endif()

  if (ShaPo_FOUND)
    include_directories( ${SHAPO_INCLUDE_DIRS} )  
    include( ${ShaPo_CONFIG} )
    include( ${SHAPO_USE_FILE} )
    list(APPEND VORO_LIBRARIES ${SHAPO_LIBRARIES} )
    add_definitions( -DHAVE_SHAPO )
  else ()
    message( FATAL_ERROR "You need libShaPo either from TPL or system to enable voro" )
  endif()
  
  message( STATUS "Voronoi with shapo enabled" )

endif()


#------------------------------------------------------------------------------#
# Enable IO
#------------------------------------------------------------------------------#

if ( VTK_FOUND OR EXODUSII_FOUND OR TECIO_FOUND ) 
  option(ENABLE_IO "Enable I/O with third party libraries." ON)
else()
  option(ENABLE_IO "Enable I/O with third party libraries." OFF)
endif()

if ( ENABLE_IO ) 
  
  if (TECIO_FOUND) 
    list( APPEND IO_LIBRARIES ${TECIO_LIBRARY} )
    include_directories( ${TECIO_INCLUDE_DIR} )
    add_definitions( -DHAVE_TECIO )
    message( STATUS "IO with tecio enabled" )
  endif()

  if(EXODUSII_FOUND)
    include_directories( ${EXODUSII_INCLUDE_DIRS} )
    list(APPEND IO_LIBRARIES ${EXODUSII_LIBRARIES} )
    add_definitions( -DHAVE_EXODUS )
    message( STATUS "IO with exodus enabled" )
  endif()

  if ( VTK_FOUND )
    if (NOT VTK_INCLUDED) 
      include( ${VTK_CONFIG} )
      include( ${VTK_USE_FILE} )
      add_definitions( -DHAVE_VTK )
    endif()
    list(APPEND IO_LIBRARIES ${VTK_LIBRARIES} )
    message( STATUS "IO with vtk enabled" )
  endif()

  if ( NOT IO_LIBRARIES )
    MESSAGE( FATAL_ERROR "You need libexodus, libtecio, or libvtk from TPL or system to enable io" )
  endif()

endif(ENABLE_IO)

#------------------------------------------------------------------------------#
# ADD to ALE libraries
#------------------------------------------------------------------------------#
list( APPEND ALE_LIBRARIES ${IO_LIBRARIES} ${VORO_LIBRARIES} )
