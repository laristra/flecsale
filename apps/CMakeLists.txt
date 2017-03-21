#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

include(regression) 

if (USE_CATALYST)
	add_subdirectory( catalyst_adaptor )
	include_directories( "${CMAKE_SOURCE_DIR}/apps/catalyst_adaptor/")
endif()

add_subdirectory( common )
add_subdirectory( hydro )
add_subdirectory( maire_hydro )

add_executable( dummy_app dummy_app.cc )
