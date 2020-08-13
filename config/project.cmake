#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# Set the minimum Cinch version
cinch_minimum_required(VERSION v1.0)

# Set the project name
project(FleCSALE)

# CMake includes
include(CMakeDependentOption)

# export compile commands to JSON file
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# We need C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)

# enable fortran
option( FLECSALE_ENABLE_FORTRAN "Enable fortran support." OFF )
if (FLECSALE_ENABLE_FORTRAN)
  enable_language(Fortran)
endif()

# Changes the Cinch defaults
option(ENABLE_BOOST_PREPROCESSOR "Enable Boost.Preprocessor" ON)

# Set header suffix regular expression
set(CINCH_HEADER_SUFFIXES "\\.h")

# set some global variables
set( FLECSALE_LIBRARIES )
set( FLECSALE_DATA_DIR "${PROJECT_SOURCE_DIR}/data" CACHE INTERNAL "")  
set( FLECSALE_TOOL_DIR "${PROJECT_SOURCE_DIR}/tools" CACHE INTERNAL "")  


#------------------------------------------------------------------------------#
# FleCSI Library
#
# Ristra libraries come first in case any options depened on what we found.
#------------------------------------------------------------------------------#

set( FLECSI_SUBMODULE_DIR ${CMAKE_SOURCE_DIR}/flecsi )
file(GLOB _flecsi_contents ${FLECSI_SUBMODULE_DIR}/*)

if ( _flecsi_contents )
  if ( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME )
    add_subdirectory( ${FLECSI_SUBMODULE_DIR} )
  endif()
  include_directories( ${FLECSI_SUBMODULE_DIR} )
  list(APPEND FLECSALE_LIBRARIES FleCSI)
  set(_runtime_path ${FleCSI_SOURCE_DIR}/flecsi/execution/${FLECSI_RUNTIME_MODEL})
  set(FLECSALE_RUNTIME_DRIVER ${_runtime_path}/runtime_driver.cc)
  set(FLECSALE_RUNTIME_MAIN   ${_runtime_path}/runtime_main.cc)
  get_property(_include_dirs DIRECTORY ${FLECSI_SUBMODULE_DIR}
    PROPERTY INCLUDE_DIRECTORIES)
  include_directories( ${_include_dirs} )
else()
  find_package(FleCSI CONFIG REQUIRED)
  include_directories(${FleCSI_INCLUDE_DIRS})
  list(APPEND FLECSALE_LIBRARIES ${FleCSI_LIBRARIES})
  set(FLECSALE_SP_RUNTIME_DRIVER ${FLECSI_RUNTIME_DRIVER})
  set(FLECSALE_SP_RUNTIME_MAIN   ${FLECSI_RUNTIME_MAIN})
endif()


set( FLECSALE_RUNTIME_MODEL ${FLECSI_RUNTIME_MODEL} )

if ( FLECSALE_RUNTIME_MODEL STREQUAL "mpi" )
  set( ENABLE_MPI ON CACHE BOOL "" FORCE)
  set( FLECSALE_UNIT_POLICY MPI )
elseif ( FLECSALE_RUNTIME_MODEL STREQUAL "legion" )
  set( FLECSALE_UNIT_POLICY LEGION )
else()
  MESSAGE( FATAL_ERROR 
    "Unknown FLECSI_SP_RUNTIME_MODEL being used: ${FLECSI_SP_RUNTIME_MODEL}" )
endif()

#------------------------------------------------------------------------------#
# Ristra Library
#------------------------------------------------------------------------------#

set( RISTRA_SUBMODULE_DIR ${CMAKE_SOURCE_DIR}/ristra )
file(GLOB _ristra_contents ${RISTRA_SUBMODULE_DIR}/*)

if (_ristra_contents )
  if ( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME )
    add_subdirectory( ${RISTRA_SUBMODULE_DIR} )
  endif()
  include_directories( ${RISTRA_SUBMODULE_DIR} )
  list(APPEND FLECSALE_LIBRARIES Ristra)
  get_property(_include_dirs DIRECTORY ${RISTRA_SUBMODULE_DIR}
    PROPERTY INCLUDE_DIRECTORIES)
  include_directories( ${_include_dirs} )
else()
  find_package(Ristra CONFIG REQUIRED)
  include_directories(${RISTRA_INCLUDE_DIRS})
  list(APPEND FLECSALE_LIBRARIES ${RISTRA_LIBRARIES})
endif()

#------------------------------------------------------------------------------#
# FleCSI-SP Library
#------------------------------------------------------------------------------#

file(GLOB _specializations_contents specializations/*)

if ( _specializations_contents )
  if ( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME )
    add_subdirectory( ${CMAKE_SOURCE_DIR}/specializations )
  endif()
  include_directories( ${CMAKE_SOURCE_DIR}/specializations )
  list(APPEND FLECSALE_LIBRARIES FleCSI-SP)
  set(FLECSI_SP_BURTON_SPECIALIZATION_INIT
    ${FleCSI-SP_SOURCE_DIR}/flecsi-sp/burton/burton_specialization_init.cc)
else()
  find_package(FleCSI-SP CONFIG REQUIRED)
  include_directories(${FleCSI_SP_INCLUDE_DIRS})
  list(APPEND FLECSALE_LIBRARIES ${FLECSI_SP_LIBRARIES})
endif()



#------------------------------------------------------------------------------#
# Begin project options
#------------------------------------------------------------------------------#

# double or single precision
set( FLECSALE_DOUBLE_PRECISION ${FLECSI_SP_DOUBLE_PRECISION} CACHE BOOL "" FORCE)

if ( FLECSALE_DOUBLE_PRECISION )
  message(STATUS "Note: Double precision build activated.")
  set( FLECSALE_TEST_TOLERANCE 1.0e-14 CACHE STRING "The testing tolerance" )
else()
  message(STATUS "Note: Sincle precision build activated.")
  set( FLECSALE_TEST_TOLERANCE 1.0e-6 CACHE STRING "The testing tolerance" )
endif()


# size of integer ids to use
set( FLECSALE_USE_64BIT_IDS ${FLECSI_SP_USE_64BIT_IDS} CACHE BOOL "" FORCE )

if( FLECSALE_USE_64BIT_IDS ) 
  message(STATUS "Note: using 64 bit integer ids.")
else()
  message(STATUS "Note: using 32 bit integer ids.")
endif()

#------------------------------------------------------------------------------#
# Enable Regression Tests
#------------------------------------------------------------------------------#

# find python for running regression tests
if(NOT TARGET Python3::Interpreter)
  find_package(Python3 COMPONENTS Interpreter)
endif()

if (Python3_Interpreter_FOUND)
  cmake_dependent_option( 
    ENABLE_REGRESSION_TESTS "Enable regression tests" ON 
    "ENABLE_UNIT_TESTS" OFF 
  )
else ()
  option(ENABLE_REGRESSION_TESTS "Enable regression tests" OFF)
endif ()

if(ENABLE_REGRESSION_TESTS AND NOT Python3_Interpreter_FOUND)
  message(FATAL_ERROR "Regression tests requested, but python was not found")
endif()

if (ENABLE_REGRESSION_TESTS)
  message (STATUS "Found PythonInterp: ${Python3_EXECUTABLE}")
endif ()

#------------------------------------------------------------------------------#
# Enable Embedded Interpreters
#------------------------------------------------------------------------------#

# find python for embedding
if(NOT TARGET Python3::Python)
  find_package(Python3 COMPONENTS Development)
endif()

option(FLECSALE_ENABLE_PYTHON "Enable Python Support" Python3_FOUND)

if (FLECSALE_ENABLE_PYTHON AND NOT Python3_FOUND)
  message(FATAL_ERROR "Python requested, but not found")
endif()

if (FLECSALE_ENABLE_PYTHON)
   message (STATUS "Found PythonLibs: ${Python3_INCLUDE_DIRS}")
   include_directories( ${Python3_INCLUDE_DIRS} )
   list( APPEND FLECSALE_LIBRARIES Python3::Python)
endif ()

# find lua for embedding
find_package (Lua 5.2 QUIET)

option(FLECSALE_ENABLE_LUA "Enable Lua Support" ${LUA_FOUND})

if (FLECSALE_ENABLE_LUA AND NOT LUA_FOUND)
  message(FATAL_ERROR "Lua requested, but not found")
endif()

if (FLECSALE_ENABLE_LUA)
   message (STATUS "Found Lua: ${LUA_INCLUDE_DIR}")
   include_directories( ${LUA_INCLUDE_DIR} )
   list( APPEND FLECSALE_LIBRARIES ${LUA_LIBRARIES} )
endif ()

#------------------------------------------------------------------------------#
# Catalyst
#------------------------------------------------------------------------------#

option(FLECSALE_ENABLE_CATALYST "Link the sim with Catalyst for in situ" OFF)

if (FLECSALE_ENABLE_CATALYST)
  find_package(ParaView REQUIRED COMPONENTS vtkPVPythonCatalyst)
  
  message(STATUS "Found Paraview: ${ParaView_DIR}")
  message(STATUS "IO with Paraview Catalyst enabled" )
  
  include("${PARAVIEW_USE_FILE}")

  if (NOT PARAVIEW_USE_MPI)
    message(SEND_ERROR "ParaView must be built with MPI enabled")
  endif()

  list( APPEND FLECSALE_LIBRARIES vtkPVPythonCatalyst vtkParallelMPI )
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

option(FLECSALE_ENABLE_VTK "Enable I/O with vtk." ${VTK_FOUND})

if(FLECSALE_ENABLE_VTK AND NOT VTK_FOUND)
  message(FATAL_ERROR "VTK requested, but not found")
endif()

if(FLECSALE_ENABLE_VTK OR FLECSALE_ENABLE_SHAPO OR FLECSALE_ENABLE_CATALYST)
    include( ${VTK_CONFIG} )
    include( ${VTK_USE_FILE} )
    list(APPEND FLECSALE_LIBRARIES ${VTK_LIBRARIES} )
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
  option(FLECSALE_ENABLE_SHAPO "Enable vornoi with shapo." ON)
else()
  option(FLECSALE_ENABLE_SHAPO "Enable vornoi with shapo." OFF)
endif()

if (FLECSALE_ENABLE_SHAPO)
  
  if (NOT VTK_FOUND)
    message( FATAL_ERROR "You need libVTK either from TPL or system to enable voro")
  endif()

  if (ShaPo_FOUND)
    include_directories( ${SHAPO_INCLUDE_DIRS} )  
    include( ${ShaPo_CONFIG} )
    include( ${SHAPO_USE_FILE} )
    list(APPEND FLECSALE_LIBRARIES ${SHAPO_LIBRARIES} )
  else ()
    message( FATAL_ERROR "You need libShaPo either from TPL or system to enable voro" )
  endif()
  
  message( STATUS "Voronoi with shapo enabled" )

endif()

#------------------------------------------------------------------------------#
# Enable Exodus
#------------------------------------------------------------------------------#

find_package(EXODUSII QUIET)

option(FLECSALE_ENABLE_EXODUS "Enable I/O with exodus." ${EXODUSII_FOUND})

if(FLECSALE_ENABLE_EXODUS AND NOT EXODUSII_FOUND)
  message(FATAL_ERROR "Exodus requested, but not found")
endif()

if(FLECSALE_ENABLE_EXODUS)
  include_directories( ${EXODUSII_INCLUDE_DIRS} )
  list(APPEND FLECSALE_LIBRARIES ${EXODUSII_LIBRARIES} )
  message( STATUS "IO with exodus enabled" )
endif()

#------------------------------------------------------------------------------#
# Boost - Right now, only used by portage
#------------------------------------------------------------------------------#

find_package(Boost 1.59.0 COMPONENTS program_options REQUIRED)

if(Boost_FOUND)
  message(STATUS "Boost includes: ${Boost_INCLUDE_DIRS}")
  message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
  include_directories( ${Boost_INCLUDE_DIRS} )
endif()

list(APPEND FLECSALE_LIBRARIES ${Boost_LIBRARIES} Boost::program_options Boost::boost)

#------------------------------------------------------------------------------#
# Portage
#------------------------------------------------------------------------------#

find_package(PORTAGE QUIET)

option(FLECSALE_ENABLE_PORTAGE "Enable Portage Support" ${PORTAGE_FOUND})

if(FLECSALE_ENABLE_PORTAGE)
  if(NOT Boost_FOUND)
    message( FATAL_ERROR "Boost is needed for Portage" )
  endif()
  message( STATUS "Portage location: ${PORTAGE_INCLUDE_DIRS}" )
  include_directories(${PORTAGE_INCLUDE_DIRS})
  list( APPEND FLECSALE_LIBRARIES ${PORTAGE_LIBRARIES} )
endif()

#------------------------------------------------------------------------------#
# Legion / MPI
#------------------------------------------------------------------------------#

find_package(Legion)

if (Legion_FOUND) 
  include_directories(${Legion_INCLUDE_DIRS})
endif()

find_package(MPI COMPONENTS C CXX REQUIRED)

if (MPI_FOUND) 
  list( APPEND FLECSALE_LIBRARIES MPI::MPI_CXX MPI::MPI_C)
endif()

#------------------------------------------------------------------------------#
# now load the extra functionality
#------------------------------------------------------------------------------#

cinch_load_extras()

#------------------------------------------------------------------------------#
# configure header
#------------------------------------------------------------------------------#

configure_file(${PROJECT_SOURCE_DIR}/config/flecsale-config.h.in
  ${CMAKE_BINARY_DIR}/flecsale-config.h @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/flecsale-config.h DESTINATION include)

include_directories(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------#
# Add library targets
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
cinch_add_library_target(FleCSALE flecsale EXPORT_TARGET FleCSALETargets)
cinch_target_link_libraries( FleCSALE ${FLECSALE_LIBRARIES} )

#------------------------------------------------------------------------------#
# Set application directory
#------------------------------------------------------------------------------#

add_subdirectory(apps)

#------------------------------------------------------------------------------#
# Extract all project options so they can be exported to the ProjectConfig.cmake
# file.
#------------------------------------------------------------------------------#

get_cmake_property(_variableNames VARIABLES)
string (REGEX MATCHALL "(^|;)FLECSALE_[A-Za-z0-9_]*" _matchedVars
  "${_variableNames}")
foreach (_variableName ${_matchedVars})
  set( FLECSALE_CONFIG_CODE
    "${FLECSALE_CONFIG_CODE}
set(${_variableName} \"${${_variableName}}\")"
  )
endforeach()

#------------------------------------------------------------------------------#
# Prepare variables for ProjectConfig file.
#------------------------------------------------------------------------------#

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY INCLUDE_DIRECTORIES)

foreach(dir ${dirs})
  if(NOT ${dir} MATCHES ${CMAKE_CURRENT_SOURCE_DIR})
    list(APPEND FLECSALE_EXTERNAL_INCLUDE_DIRS ${dir})
  endif()
endforeach()

set(FLECSALE_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR})
set(FLECSALE_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(FLECSALE_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSALE)

#------------------------------------------------------------------------------#
# Export targets and package.
#------------------------------------------------------------------------------#

export(
  TARGETS FleCSALE
  FILE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSALETargets.cmake
)

export(PACKAGE FleCSALE)

#------------------------------------------------------------------------------#
# configure .cmake file (for other projects)
#------------------------------------------------------------------------------#

configure_file(${PROJECT_SOURCE_DIR}/config/FleCSALEConfig.cmake.in
  ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSALEConfig.cmake @ONLY)

install(
  FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSALEConfig.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSALE
)

install(
  EXPORT FleCSALETargets
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSALE
  COMPONENT dev
)

