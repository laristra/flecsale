#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro main.cc )
target_link_libraries( hydro ${ALE_LIBRARIES} )


create_regression_test( 
	NAME shock_box 
	COMMAND $<TARGET_FILE:hydro> 
	COMPARE shock_box_3d0000048.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000048.dat.std 
)
