#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(io_HEADERS
  catalyst/adaptor.h
  io_exodus.h
  write_binary.h
  vtk.h
  
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

set(io_SOURCES
  vtk.cc
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

