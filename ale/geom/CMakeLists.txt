#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(geom_HEADERS
  centroid.h    detail/centroid.h
  area.h        detail/area.h
  geometric_shapes.h
  normal.h
  point.h
  PARENT_SCOPE
)

cinch_add_unit(geom
    SOURCES 
      test/centroid.cc
)


#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
