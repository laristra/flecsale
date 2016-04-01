#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set( ALE_LIBS ale flecsi ${IO_LIBRARIES} ${SHAPO_LIBRARIES} )

add_subdirectory( hydro )
add_subdirectory( lag_hydro )

add_executable( dummy_app dummy_app.cc )


#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
