#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro main.cc )
target_link_libraries( hydro ${ALE_LIBRARIES} )
