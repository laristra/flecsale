#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(geom_HEADERS
  centroid.h    detail/centroid_impl.h
  area.h        detail/area_impl.h
  normal.h
  point.h
  shapes/hexahedron.h
  shapes/quadrilateral.h
  shapes/polyhedron.h
  shapes/polygon.h
  shapes/tetrahedron.h
  shapes/triangle.h
  PARENT_SCOPE
)

cinch_add_unit(geom
    SOURCES 
      test/centroid.cc
      shapes/test/shapes.cc
)

