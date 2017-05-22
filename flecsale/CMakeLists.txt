#------------------------------------------------------------------------------#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# create a list to append object libraries to
set( FleCSALE_OBJECTS )

# add the individual submodules
add_subdirectory( common )
add_subdirectory( eos )
add_subdirectory( eqns )
add_subdirectory( geom )
add_subdirectory( io )
add_subdirectory( linalg )
add_subdirectory( math )
add_subdirectory( mesh )
add_subdirectory( utils )
add_subdirectory( fortran )

# in order to enable shared libraries, they need to be built using fPIC
set_target_properties( 
  ${FleCSALE_OBJECTS} PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)

# turn the list of target names into a list of objects
set(objects)
foreach( obj ${FleCSALE_OBJECTS} )
  list(APPEND objects $<TARGET_OBJECTS:${obj}>)
endforeach()


# create the library and link it to thirparties
add_library(flecsale ${objects})
target_link_libraries( flecsale ${FleCSALE_LIBRARIES} )

# install the library
install( TARGETS flecsale DESTINATION lib)

#~---------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#~---------------------------------------------------------------------------~-#
