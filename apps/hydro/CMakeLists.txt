#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro main.cc )
target_link_libraries( hydro ${ALE_LIBRARIES} )


create_regression_test( shock_box hydro ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000033.dat.std shock_box_3d0000033.dat )
