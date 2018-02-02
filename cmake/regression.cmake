#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

# required includes
include(ProcessorCount)

# get the scripts directory
set(REGRESSION_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR} )

#-------------------------------------------------------------------------------
# This macro creates a regression test

function(create_regression_test)
  if (ENABLE_REGRESSION_TESTS)

    # the command to run to compare outputs
    set (TEST_COMMAND "${PYTHON_EXECUTABLE} ${FLECSALE_TOOL_DIR}/numdiff.py --verbose --absolute ${TEST_TOLERANCE}")

    # parse the arguments
    set(options)
    set(oneValueArgs NAME COMPARE STANDARD THREADS)
    set(multiValueArgs COMMAND INPUTS)
    cmake_parse_arguments(args "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
    # check the preconditions
    if( NOT args_NAME )
      message( FATAL_ERROR "You must specify a test name using NAME." )
    endif()
  
    if( NOT args_COMMAND )
      message( FATAL_ERROR "You must specify a test command using COMMAND." )
    endif()
  
    if( NOT args_COMPARE )
      message( FATAL_ERROR "You must specify a file to compare against the "
        "standard with using COMPARE.")
    endif()
  
    if( NOT args_COMPARE )
      message( FATAL_ERROR "You must specify a standard file for comparisons "
        "using STANDARD.")
    endif()
  
    # use at least one thread
    if( NOT args_THREADS )
      set( args_THREADS 1 )
    endif()
  
    # count number of processors
    set( NUM_PROCS 1 )
    if ( ENABLE_OPENMP )
      set( NUM_PROCS OMP_NUM_PROCS )
    endif()    
  
    if ( NOT( NUM_PROCS LESS args_THREADS) )
 
      # add the test
      add_test( 
        NAME ${args_NAME}
        COMMAND ${CMAKE_COMMAND}
          "-Dtest_cmd=${args_COMMAND}"
          -Dcompare_cmd=${TEST_COMMAND}
          -Doutput_blessed=${args_STANDARD}
          -Doutput_test=${args_COMPARE}
          -P ${REGRESSION_CMAKE_DIR}/run_test.cmake
      )
  
      # for openmp
      SET_TESTS_PROPERTIES( ${args_NAME}
        PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=${args_THREADS}")
      
    endif()
  
  endif(ENABLE_REGRESSION_TESTS)
endfunction()
 
