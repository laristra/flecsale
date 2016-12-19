#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

include(regression) 

add_definitions( -DFLECSI_RUNTIME_MODEL_serial )

add_subdirectory( hydro )
add_subdirectory( maire_hydro )

add_executable( dummy_app dummy_app.cc )
