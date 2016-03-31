#~----------------------------------------------------------------------------~#
#  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
# /@@/////  /@@          @@////@@ @@////// /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@
# /@@       /@@/@@//// //@@    @@       /@@/@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
# //       ///  //////   //////  ////////  // 
# 
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(mesh_HEADERS
  burton/burton_entity_types.h
  burton/burton.h
  burton/burton_io_exodus.h
  burton/burton_io_tecplot.h
  burton/burton_io_vtk.h
  burton/burton_mesh.h
  burton/burton_mesh_traits.h
  burton/burton_types.h
  factory.h
  voronoi.h  detail/voronoi.h
  PARENT_SCOPE
)


set(mesh_SOURCES
  burton/burton_entity_types.cc
  PARENT_SCOPE
)

set( mesh_LIBS 
  ${IO_LIBRARIES} 
  ${SHAPO_LIBRARIES}
  flecsi 
)

cinch_add_unit(mesh
    SOURCES 
      burton/test/burton.cc
      burton/test/burton_copy.cc
      burton/test/burton_io.cc
      burton/test/gradients.cc
      burton/test/voronoi.cc
    INPUTS 
      burton/test/burton.blessed
      burton/test/mixed/mixed.g
    LIBRARIES ${mesh_LIBS} ale
)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
