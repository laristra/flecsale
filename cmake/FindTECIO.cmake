#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find libtecio
# Find the native TECIO headers and libraries.
#
#  TECIO_INCLUDE_DIRS - where to find exodusII.h, etc.
#  TECIO_LIBRARIES    - List of libraries when using xoIIv2c.
#  TECIO_FOUND        - True if exodus found.
#

find_path(TECIO_INCLUDE_DIR TECIO.h)

find_library(TECIO_LIBRARY NAMES tecio)

set(TECIO_LIBRARIES ${TECIO_LIBRARY} )
set(TECIO_INCLUDE_DIRS ${TECIO_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TECIO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TECIO DEFAULT_MSG TECIO_LIBRARY TECIO_INCLUDE_DIR )

mark_as_advanced(TECIO_INCLUDE_DIR TECIO_LIBRARY)
