#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 2.8)

include(CMakeDependentOption)

#------------------------------------------------------------------------------#
# Set the project name
#------------------------------------------------------------------------------#

project(FleCSALE)

# enable fortran
option( ENABLE_FORTRAN "Enable fortran support." OFF )
if (ENABLE_FORTRAN)
  enable_language(Fortran)
endif()

#------------------------------------------------------------------------------#
# Include project-level CMake configuration file
#------------------------------------------------------------------------------#

if(NOT "$ENV{CINCH_SOURCE_DIR}" STREQUAL "")
  set(CINCH_SOURCE_DIR $ENV{CINCH_SOURCE_DIR})
elseif(EXISTS ${CMAKE_SOURCE_DIR}/cinch)
  set(CINCH_SOURCE_DIR ${CMAKE_SOURCE_DIR}/cinch)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CINCH_SOURCE_DIR}/cmake")
include(common_includes)
include(install_headers)

# This changes the Cinch default
option(ENABLE_BOOST_PREPROCESSOR "Enable Boost.Preprocessor" ON)

# now load the extra functionality
cinch_load_extras()

#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#

cinch_minimum_required(2.0)

#------------------------------------------------------------------------------#
# Set the C++ standard
#------------------------------------------------------------------------------#

# We need C++ 14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)

#------------------------------------------------------------------------------#
# Begin project setup
#------------------------------------------------------------------------------#

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# set some global variables
set( FleCSALE_LIBRARIES )
set( FleCSALE_DATA_DIR "${PROJECT_SOURCE_DIR}/data" CACHE INTERNAL "")  
set( FleCSALE_TOOL_DIR "${PROJECT_SOURCE_DIR}/tools" CACHE INTERNAL "")  

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
   list( APPEND FleCSALE_LIBRARIES ${PYTHON_LIBRARIES} )
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
   list( APPEND FleCSALE_LIBRARIES ${LUA_LIBRARIES} )
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
  list( APPEND FleCSALE_LIBRARIES ${OPENSSL_LIBRARIES} )
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
  list( APPEND FleCSALE_LIBRARIES ${Caliper_LIBRARIES} )
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

  list( APPEND FleCSALE_LIBRARIES vtkPVPythonCatalyst vtkParallelMPI )
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
    #include( ${VTK_CONFIG} )
    include( ${VTK_USE_FILE} )
    add_definitions( -DHAVE_VTK )
    list(APPEND FleCSALE_LIBRARIES ${VTK_LIBRARIES} )
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
    list(APPEND FleCSALE_LIBRARIES ${SHAPO_LIBRARIES} )
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
  list( APPEND FleCSALE_LIBRARIES ${TECIO_LIBRARY} )
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
  list(APPEND FleCSALE_LIBRARIES ${EXODUSII_LIBRARIES} )
  message( STATUS "IO with exodus enabled" )
endif()

#------------------------------------------------------------------------------#
# Boost - Right now, only used by portage
#------------------------------------------------------------------------------#

find_package(Boost QUIET)

if(Boost_FOUND)
  message(STATUS "Boost location: ${Boost_INCLUDE_DIRS}")
  include_directories( ${Boost_INCLUDE_DIRS} )
endif()

#------------------------------------------------------------------------------#
# Portage
#------------------------------------------------------------------------------#

find_package(PORTAGE QUIET)

option(ENABLE_PORTAGE "Enable Portage Support" ${PORTAGE_FOUND})

if(ENABLE_PORTAGE)
  if(NOT Boost_FOUND)
    message( FATAL_ERROR "Boost is needed for Portage" )
  endif()

  find_package(LAPACKE REQUIRED)
  include_directories(${LAPACKE_INCLUDE_DIRS})
  add_definitions("-DHAVE_LAPACKE")

  include_directories(${PORTAGE_INCLUDE_DIR})
  add_definitions(-DHAVE_PORTAGE -DPORTAGE_SERIAL_ONLY)
  list( APPEND FleCSALE_LIBRARIES ${PORTAGE_LIBRARY} ${LAPACKE_LIBRARIES} )
endif()


#------------------------------------------------------------------------------#
# Add subprojects
#------------------------------------------------------------------------------#

set(FLECSI_RUNTIME_MODEL "serial" CACHE STRING
  "Select the runtime model")
set_property(CACHE FLECSI_RUNTIME_MODEL
  PROPERTY STRINGS serial mpilegion legion mpi)

find_package(FleCSI QUIET)
if(NOT FleCSI_FOUND)
  add_subdirectory( flecsi )
  set( FleCSI_LIBRARIES flecsi )
  include_directories( flecsi )
  include_directories( flecsi/flecsi )
  include_directories( ${CMAKE_BINARY_DIR} )
else()
  include_directories(${FleCSI_INCLUDE_DIRS})
  list( APPEND FleCSALE_LIBRARIES ${FleCSI_LIBRARIES})
endif()
list(APPEND FleCSALE_LIBRARIES ${FleCSI_LIBRARIES} )

set(FleCSI_EXEC_DIR ${PROJECT_SOURCE_DIR}/flecsi/flecsi/execution)
set(
  FleCSI_RUNTIME_MAIN ${FleCSI_EXEC_DIR}/runtime_main.cc 
  CACHE INTERNAL ""
)

set(FleCSI_RUNTIME_DRIVER 
  ${FleCSI_EXEC_DIR}/${FLECSI_RUNTIME_MODEL}/runtime_driver.cc
  CACHE INTERNAL ""
)
add_definitions( -DFLECSI_RUNTIME_MODEL_${FLECSI_RUNTIME_MODEL} )

#------------------------------------------------------------------------------#
# Add library targets
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(flecsale)

#------------------------------------------------------------------------------#
# Set application directory
#------------------------------------------------------------------------------#

add_subdirectory(apps)

