#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro_2d main.cc inputs.cc )
target_link_libraries( hydro_2d common ${ALE_LIBRARIES} )


create_regression_test( 
	NAME shock_box_2d
	COMMAND $<TARGET_FILE:hydro_2d> -f ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_2d.lua
	COMPARE shock_box_2d0000007.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_2d0000007.dat.std 
)

create_regression_test( 
	NAME shock_box_2d_omp4
	COMMAND $<TARGET_FILE:hydro_2d> -f ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_2d.lua
	THREADS 4
	COMPARE shock_box_2d0000007.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_2d0000007.dat.std 
)
