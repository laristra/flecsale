#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set( LIBS ale flecsi ${IO_LIBRARIES} )

add_executable( hydro hydro.cc )
target_link_libraries( hydro ${LIBS} )

add_executable( dummy_app dummy_app.cc )


#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
