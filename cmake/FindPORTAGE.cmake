#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find portage
# Find the native Portage headers and libraries.
#
#  PORTAGE_INCLUDE_DIRS - where to find flecsi.h, etc.
#  PORTAGE_LIBRARIES    - List of libraries when using flecsi.
#  PORTAGE_FOUND        - True if flecsi found.

# Look for the header file.
FIND_PATH(PORTAGE_INCLUDE_DIR NAMES portage.h
  PATH_SUFFIXES portage/support)

if ( PORTAGE_INCLUDE_DIR )
  set(PORTAGE_INCLUDE_DIR ${PORTAGE_INCLUDE_DIR}/../.. )
endif()

# Look for the library.
FIND_LIBRARY(PORTAGE_LIBRARY NAMES portage libportage)
FIND_LIBRARY(WONTON_LIBRARY NAMES wonton libwonton)

# handle the QUIETLY and REQUIRED arguments and set PORTAGE_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PORTAGE PORTAGE_LIBRARY PORTAGE_INCLUDE_DIR)

# Copy the results to the output variables.
SET(PORTAGE_LIBRARIES ${PORTAGE_LIBRARY} ${WONTON_LIBRARY})
SET(PORTAGE_INCLUDE_DIRS ${PORTAGE_INCLUDE_DIR})

MARK_AS_ADVANCED(PORTAGE_INCLUDE_DIR PORTAGE_LIBRARY)
