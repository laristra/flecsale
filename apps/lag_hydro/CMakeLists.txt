#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( lag_hydro main.cc )
target_link_libraries( lag_hydro ${ALE_LIBRARIES} )

create_regression_test( 
	NAME sedov
	COMMAND $<TARGET_FILE:lag_hydro>
	COMPARE sedov_3d0000010.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
)

create_regression_test( 
	NAME sedov_omp4
	COMMAND $<TARGET_FILE:lag_hydro>
	THREADS 4
	COMPARE sedov_3d0000010.dat 
	STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
)

