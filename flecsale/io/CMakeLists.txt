#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(io_HEADERS
  catalyst/adaptor.h
  write_binary.h
  vtk.h
)

set(io_SOURCES
  vtk.cc
)

mcinch_install_headers(
  FILES ${io_HEADERS}
  DESTINATION include/flecsale/io
)

add_library(flecsale_io OBJECT ${io_SOURCES})
set(FleCSALE_OBJECTS ${FleCSALE_OBJECTS} flecsale_io PARENT_SCOPE)
