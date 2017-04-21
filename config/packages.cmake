#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

include(CMakeDependentOption)

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
# Catalyst
#------------------------------------------------------------------------------#

option(ENABLE_CATALYST "Link the sim with Catalyst for in situ" OFF)

if (ENABLE_CATALYST)
  find_package(ParaView REQUIRED COMPONENTS vtkPVPythonCatalyst)
  
  message(STATUS "Found Paraview: ${ParaView_DIR}")
  message(STATUS "IO with Paraview Catalyst enabled" )
  
  include("${PARAVIEW_USE_FILE}")
  add_definitions(-DHAVE_CATALYST)
  #add_definitions( -DHAVE_VTK )

  if (NOT PARAVIEW_USE_MPI)
    message(SEND_ERROR "ParaView must be built with MPI enabled")
  endif()

  list( APPEND ALE_LIBRARIES vtkPVPythonCatalyst vtkParallelMPI )
endif()


#------------------------------------------------------------------------------#
# Enable VTK
#------------------------------------------------------------------------------#

# if VTK still has not been found by another package, try to find it directly
if (NOT VTK_FOUND)
  find_package(VTK QUIET)
  if (VTK_FOUND)
    message(STATUS "Found VTK: ${VTK_DIR}")
  endif()
endif()

option(ENABLE_VTK "Enable I/O with vtk." ${VTK_FOUND})

if(ENABLE_VTK AND NOT VTK_FOUND)
  message(FATAL_ERROR "VTK requested, but not found")
endif()

if(ENABLE_VTK OR ENABLE_SHAPO OR ENABLE_CATALYST)
    include( ${VTK_CONFIG} )
    include( ${VTK_USE_FILE} )
    add_definitions( -DHAVE_VTK )
    list(APPEND ALE_LIBRARIES ${VTK_LIBRARIES} )
    message( STATUS "IO with vtk enabled" )
endif()


#------------------------------------------------------------------------------#
# Enable shapo
#------------------------------------------------------------------------------#

find_package(ShaPo QUIET)
if (ShaPo_FOUND)
  message(STATUS "Found ShaPo: ${ShaPo_DIR}")
endif()


if ( ShaPo_FOUND AND VTK_FOUND )
  option(ENABLE_SHAPO "Enable vornoi with shapo." ON)
else()
  option(ENABLE_SHAPO "Enable vornoi with shapo." OFF)
endif()

if(ENABLE_SHAPO)
  
  if (NOT VTK_FOUND)
    message( FATAL_ERROR "You need libVTK either from TPL or system to enable voro")
  endif()

  if (ShaPo_FOUND)
    include_directories( ${SHAPO_INCLUDE_DIRS} )  
    include( ${ShaPo_CONFIG} )
    include( ${SHAPO_USE_FILE} )
    list(APPEND ALE_LIBRARIES ${SHAPO_LIBRARIES} )
    add_definitions( -DHAVE_SHAPO )
  else ()
    message( FATAL_ERROR "You need libShaPo either from TPL or system to enable voro" )
  endif()
  
  message( STATUS "Voronoi with shapo enabled" )

endif()

#------------------------------------------------------------------------------#
# Enable Tecio
#------------------------------------------------------------------------------#

find_package(TECIO QUIET)


if(ENABLE_TECIO AND NOT TECIO_FOUND)
  message(FATAL_ERROR "Tecplot requested, but not found")
endif()

if(ENABLE_TECIO)
  include_directories( ${TECIO_INCLUDE_DIR} )
  add_definitions( -DHAVE_TECIO )
  list( APPEND ALE_LIBRARIES ${TECIO_LIBRARY} )
  message( STATUS "IO with tecio enabled" )
endif()

#------------------------------------------------------------------------------#
# Enable Exodus
#------------------------------------------------------------------------------#

find_package(EXODUSII QUIET)

option(ENABLE_EXODUS "Enable I/O with exodus." ${EXODUSII_FOUND})

if(ENABLE_EXODUS AND NOT EXODUSII_FOUND)
  message(FATAL_ERROR "Exodus requested, but not found")
endif()

if(ENABLE_EXODUS)
  include_directories( ${EXODUSII_INCLUDE_DIRS} )
  add_definitions( -DHAVE_EXODUS )
  list(APPEND ALE_LIBRARIES ${EXODUSII_LIBRARIES} )
  message( STATUS "IO with exodus enabled" )
endif()

