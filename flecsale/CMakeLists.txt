#------------------------------------------------------------------------------#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Add public library headers
#------------------------------------------------------------------------------#

add_subdirectory( common )
add_subdirectory( eos )
add_subdirectory( eqns )
add_subdirectory( geom )
add_subdirectory( io )
add_subdirectory( linalg )
add_subdirectory( math )
add_subdirectory( mesh )
add_subdirectory( utils )

set(objects
  flecsale_common
  flecsale_eos
  flecsale_eqns
  flecsale_geom
  flecsale_io
  flecsale_linalg
  flecsale_math
  flecsale_mesh
  flecsale_utils
)

set(valid_objects)
foreach( obj ${objects} )
  if ( TARGET ${obj} )
    set(valid_objects ${valid_objects} $<TARGET_OBJECTS:${obj}>)
  endif()
endforeach()

add_library(flecsale ${valid_objects})
target_link_libraries( flecsale ${FleCSALE_LIBRARIES} )
install( TARGETS flecsale DESTINATION lib)

#~---------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#~---------------------------------------------------------------------------~-#
