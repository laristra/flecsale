#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( maire_hydro main.cc )
target_link_libraries( maire_hydro ${ALE_LIBRARIES} )

create_regression_test( 
	NAME sedov_maire
	COMMAND $<TARGET_FILE:maire_hydro> --postfix dat --max_steps 10 
	  --output_freq 10
	COMPARE sedov_3d0000010.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
)

create_regression_test( 
	NAME sedov_maire_omp4
	COMMAND $<TARGET_FILE:maire_hydro> --postfix dat --max_steps 10
	  --output_freq 10
	THREADS 4
	COMPARE sedov_3d0000010.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
)

