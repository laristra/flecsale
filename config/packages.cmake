#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

include(CMakeDependentOption)

#------------------------------------------------------------------------------#
# Global variables
#------------------------------------------------------------------------------#

set( IO_LIBRARIES )
set( VORO_LIBRARIES )

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

# flecsi needs Cereal
find_package( Cereal REQUIRED )
include_directories(${Cereal_INCLUDE_DIR})

#------------------------------------------------------------------------------#
# Enable Regression Tests
#------------------------------------------------------------------------------#

# find python for running regression tests
find_package (PythonInterp QUIET)
if (PYTHONINTERP_FOUND)
  cmake_dependent_option( 
    ENABLE_REGRESSION_TESTS "Enable regression tests" ON 
    "ENABLE_UNIT_TESTS" OFF 
  )
else ()
  option(ENABLE_REGRESSION_TESTS "Enable regression tests" OFF)
endif ()

if(ENABLE_REGRESSION_TESTS AND NOT PYTHONINTERP_FOUND)
  message(FATAL_ERROR "Regression tests requested, but python was not found")
endif()

if (ENABLE_REGRESSION_TESTS)
  message (STATUS "Found PythonInterp: ${PYTHON_EXECUTABLE}")
endif ()

#------------------------------------------------------------------------------#
# Enable Embedded Interpreters
#------------------------------------------------------------------------------#

# find python for embedding
find_package (PythonLibs QUIET)

option(ENABLE_PYTHON "Enable Python Support" ${PYTHONLIBS_FOUND})

if(ENABLE_PYTHON AND NOT PYTHONLIBS_FOUND)
  message(FATAL_ERROR "Python requested, but not found")
endif()

if (ENABLE_PYTHON)
   message (STATUS "Found PythonLibs: ${PYTHON_INCLUDE_DIRS}")
   include_directories( ${PYTHON_INCLUDE_DIRS} )
   list( APPEND ALE_LIBRARIES ${PYTHON_LIBRARIES} )
   add_definitions( -DHAVE_PYTHON )    
endif ()

# find lua for embedding
find_package (Lua 5.2 QUIET)

option(ENABLE_LUA "Enable Lua Support" ${LUA_FOUND})

if(ENABLE_LUA AND NOT LUA_FOUND)
  message(FATAL_ERROR "Lua requested, but not found")
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

option(ENABLE_OPENSSL "Enable OpenSSL Support" ${OPENSSL_FOUND})

if(ENABLE_OPENSSL AND NOT OPENSSL_FOUND)
  message(FATAL_ERROR "OpenSSL requested, but not found")
endif()

if(ENABLE_OPENSSL)
  message(STATUS "Found OpenSSL: ${OPENSSL_INCLUDE_DIR}")
  include_directories(${OPENSSL_INCLUDE_DIR})
  # flecsi uses ENABLE_OPENSSL, flecsale uses HAVE_OPENSSL
  add_definitions(-DHAVE_OPENSSL -DENABLE_OPENSSL)
  list( APPEND ALE_LIBRARIES ${OPENSSL_LIBRARIES} )
endif()

#------------------------------------------------------------------------------#
# Caliper
#------------------------------------------------------------------------------#

find_package(Caliper QUIET)

option(ENABLE_CALIPER "Enable Caliper Support" ${Caliper_FOUND})

if(ENABLE_CALIPER AND NOT Caliper_FOUND)
  message(FATAL_ERROR "Caliper requested, but not found")
endif()
  
if(ENABLE_CALIPER)
  message(STATUS "Found Caliper: ${Caliper_INCLUDE_DIRS}")
  include_directories(${Caliper_INCLUDE_DIRS})
  add_definitions(-DHAVE_CALIPER)
  list( APPEND ALE_LIBRARIES ${Caliper_LIBRARIES} )
endif()


#------------------------------------------------------------------------------#
# Find some general libraries
#------------------------------------------------------------------------------#

find_package(EXODUSII)

find_package(TECIO QUIET)

find_package(ShaPo QUIET)
if (ShaPo_FOUND)
  message(STATUS "Found ShaPo: ${ShaPo_DIR}")
endif()

if (NOT VTK_FOUND)
  find_package(VTK QUIET)
  if (VTK_FOUND)
    message(STATUS "Found VTK: ${VTK_DIR}")
  endif()
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


