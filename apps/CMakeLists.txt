#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#



# the command to run to compare outputs
set (TEST_COMMAND "${PYTHON_EXECUTABLE} ${ALE_TOOL_DIR}/numdiff.py --absolute ${TEST_TOLERANCE}")


#-------------------------------------------------------------------------------
# some macros for creating tests

function(create_regression_test name target std compare)

  
   cmake_parse_arguments(arg "" "" "ARGS" ${ARGN})

   # add the test
   add_test( 
             NAME ${name}
             COMMAND ${CMAKE_COMMAND}
              -Dtest_cmd=$<TARGET_FILE:${target}>
              -Dtest_args="${args_ARGS}"
              -Dcompare_cmd=${TEST_COMMAND}
              -Doutput_blessed=${std}
              -Doutput_test=${compare}
              -P ${CMAKE_SOURCE_DIR}/cmake/run_test.cmake
   )

endfunction()
 

#-------------------------------------------------------------------------------


add_subdirectory( hydro )
add_subdirectory( lag_hydro )

add_executable( dummy_app dummy_app.cc )
