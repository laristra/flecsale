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
  burton/burton_io_vtm.h
  burton/burton_io_vtu.h
  burton/burton_mesh.h
  burton/burton_mesh_traits.h
  burton/burton_types.h
  burton/burton_vtk_utils.h
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

# set some global paths
set( ALE_DATA_DIR ../../data )

cinch_add_unit(mesh
    SOURCES 
      burton/test/burton.cc
      burton/test/burton_copy.cc
      burton/test/burton_io.cc
      burton/test/burton_voro.cc
      burton/test/gradients.cc
    INPUTS 
      burton/test/burton.blessed
      ${ALE_DATA_DIR}/meshes/mixed/mixed.g
      ${ALE_DATA_DIR}/meshes/mixed/mixed.exo
      ${ALE_DATA_DIR}/meshes/mixed/mixed.vtu
      ${ALE_DATA_DIR}/meshes/mixed/mixed.vtm
      ${ALE_DATA_DIR}/meshes/mixed/mixed_0_0.vtu
      ${ALE_DATA_DIR}/meshes/mixed/mixed_1_0.vtu
      ${ALE_DATA_DIR}/meshes/mixed/mixed.vtk
      ${ALE_DATA_DIR}/meshes/voro/voro.g
    LIBRARIES ${mesh_LIBS} ale
)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
