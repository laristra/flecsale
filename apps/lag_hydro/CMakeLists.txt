#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( lag_hydro main.cc )
target_link_libraries( lag_hydro ${ALE_LIBRARIES} )
