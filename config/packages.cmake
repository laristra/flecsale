#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Flecsi Config
#------------------------------------------------------------------------------#

set(FLECSI_EXEC_DIR ${PROJECT_SOURCE_DIR}/flecsi/flecsi/execution)
set(FLECSI_RUNTIME_MAIN ${FLECSI_EXEC_DIR}/runtime_main.cc)

if(FLECSI_RUNTIME_MODEL STREQUAL "serial")
  set(FLECSI_RUNTIME_DRIVER ${FLECSI_EXEC_DIR}/serial/runtime_driver.cc)
  add_definitions( -DFLECSI_RUNTIME_MODEL_serial )
elseif(FLECSI_RUNTIME_MODEL STREQUAL "legion")
  set(FLECSI_RUNTIME_DRIVER ${FLECSI_EXEC_DIR}/legion/runtime_driver.cc)
  add_definitions( -DFLECSI_RUNTIME_MODEL_legion )
elseif(FLECSI_RUNTIME_MODEL STREQUAL "mpilegion")
  set(FLECSI_RUNTIME_DRIVER ${FLECSI_EXEC_DIR}/mpilegion/runtime_driver.cc)
  add_definitions( -DFLECSI_RUNTIME_MODEL_mpilegion )
else()
    message(FATAL_ERROR "This runtime is not yet supported")
endif()

#------------------------------------------------------------------------------#
# Global variables
#------------------------------------------------------------------------------#

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
find_package (PythonInterp)
if (PYTHONINTERP_FOUND)
   message (STATUS "Found PythonInterp: ${PYTHON_EXECUTABLE}")
else ()
  message (FATAL_ERROR "Did not find python. Python is needed to run regression tests.")
endif ()

#------------------------------------------------------------------------------#
# Enable Embedded Interpreters
#------------------------------------------------------------------------------#

# find python for embedding
find_package (PythonLibs QUIET)

if (PYTHONLIBS_FOUND)
  option(ENABLE_PYTHON "Enable Python Support" ON)
else()
  option(ENABLE_PYTHON "Enable Python Support" OFF)
endif()

if (ENABLE_PYTHON)
   message (STATUS "Found PythonLibs: ${PYTHON_INCLUDE_DIRS}")
   include_directories( ${PYTHON_INCLUDE_DIRS} )
   list( APPEND ALE_LIBRARIES ${PYTHON_LIBRARIES} )
   add_definitions( -DHAVE_PYTHON )    
endif ()

# find lua for embedding
find_package (Lua 5.2 QUIET)

if (LUA_FOUND)
  option(ENABLE_LUA "Enable Lua Support" ON)
else()
  option(ENABLE_LUA "Enable Lua Support" OFF)
endif()

if (ENABLE_LUA)
   message (STATUS "Found Lua: ${LUA_INCLUDE_DIR}")
   include_directories( ${LUA_INCLUDE_DIR} )
   list( APPEND ALE_LIBRARIES ${LUA_LIBRARIES} )
   add_definitions( -DHAVE_LUA )    
endif ()

#------------------------------------------------------------------------------#
# OpenSSL
#------------------------------------------------------------------------------#

# OpenSSL
find_package(OpenSSL QUIET)

if (OPENSSL_FOUND)
  option(ENABLE_OPENSSL "Enable OpenSSL Support" ON)
else()
  option(ENABLE_OPENSSL "Enable OpenSSL Support" OFF)
endif()

if(ENABLE_OPENSSL)
  include_directories(${OPENSSL_INCLUDE_DIR})
  add_definitions(-DHAVE_OPENSSL)
  list( APPEND ALE_LIBRARIES ${OPENSSL_LIBRARIES} )
endif()


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
list( APPEND ALE_LIBRARIES 
  ${IO_LIBRARIES} 
  ${VORO_LIBRARIES}
)


