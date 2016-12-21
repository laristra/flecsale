# some argument checking:

# test_cmd is the command to run with all its arguments
if( NOT test_cmd )
   message( FATAL_ERROR "Variable test_cmd not defined" )
endif()

# compare_cmd is the command to compare files
if( NOT compare_cmd )
   message( FATAL_ERROR "Variable compare_cmd not defined" )
endif()

# output_blessed contains the name of the "blessed" output file
if( NOT output_blessed )
   message( FATAL_ERROR "Variable output_blessed not defined" )
endif()

# output_test contains the name of the output file the test_cmd will produce
if( NOT output_test )
   message( FATAL_ERROR "Variable output_test not defined" )
endif()

# how many threads are there
message(STATUS "Using $ENV{OMP_NUM_THREADS} threads")

# blow away the compare-to-file in case it is already there
file(REMOVE ${output_test})

# run the test
separate_arguments( test_cmd ) 

string(REPLACE ";" " " test_cmd_string "${test_cmd}")
message(STATUS "Executing '${test_cmd_string}'")

execute_process(
   COMMAND ${test_cmd}
   #OUTPUT_FILE log
   RESULT_VARIABLE test_not_successful
)

if( test_not_successful )
   message( SEND_ERROR "Error running ${test_cmd}" )
endif()

# need to fix the spaces in the passed command for some reason
separate_arguments( compare_cmd ) 

# run the diff
string(REPLACE ";" " " test_cmd_string "${compare_cmd} ${output_blessed} ${output_test}")
message(STATUS "Executing '${test_cmd_string}'")

execute_process(
  COMMAND ${compare_cmd} ${output_blessed} ${output_test}
  RESULT_VARIABLE test_not_successful
)

if( test_not_successful )
   message( SEND_ERROR "${output_test} does not match ${output_blessed}!" )
endif()
