#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro main.cc inputs.cc )
target_link_libraries( hydro common ${ALE_LIBRARIES} )


create_regression_test( 
	NAME shock_box 
	COMMAND $<TARGET_FILE:hydro> --postfix dat --output_freq 11
	COMPARE shock_box_3d0000011.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000011.dat.std 
)

create_regression_test( 
	NAME shock_box_omp4
	COMMAND $<TARGET_FILE:hydro> --postfix dat --output_freq 11
	THREADS 4
	COMPARE shock_box_3d0000011.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000011.dat.std 
)
