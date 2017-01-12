#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( maire_hydro_2d main.cc inputs.cc )
target_link_libraries( maire_hydro_2d common ${ALE_LIBRARIES} )

create_regression_test( 
	NAME sedov_maire_2d
	COMMAND $<TARGET_FILE:maire_hydro_2d> -f ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d.lua
	COMPARE sedov_2d0000010.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d0000010.dat.std 
)

create_regression_test( 
	NAME sedov_maire_2d_omp4
	COMMAND $<TARGET_FILE:maire_hydro_2d> -f ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d.lua
	THREADS 4
	COMPARE sedov_2d0000010.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d0000010.dat.std 
)

