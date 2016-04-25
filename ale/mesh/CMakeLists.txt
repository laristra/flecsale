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

  burton/burton_triangle.h
  burton/burton_quadrilateral.h
  burton/burton_polygon.h

  burton/burton_tetraahedron.h
  burton/burton_hexahedron.h
  burton/burton_polyhedron.h

  factory.h
  voronoi.h  detail/voronoi.h
  vtk_utils.h

  PARENT_SCOPE
)


# set(mesh_SOURCES
#   burton/burton_entity_types.cc
#   PARENT_SCOPE
# )

set( mesh_LIBS 
  ${IO_LIBRARIES} 
  ${SHAPO_LIBRARIES}
  flecsi 
)

# set some global paths
set( ALE_DATA_DIR ../../data )

cinch_add_unit(mesh
    SOURCES 
      burton/test/burton_test_base.h
      burton/test/burton_create_test.h
      burton/test/burton_nd_test.h
      burton/test/burton_io_test.h
      burton/test/burton_voro_test.h

      burton/test/burton_create.cc
      burton/test/burton_2d.cc
      burton/test/burton_3d.cc
      burton/test/burton_io.cc
      burton/test/burton_voro.cc
      burton/test/burton_grad.cc
    INPUTS 
      burton/test/burton.blessed
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed.g
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed.vtu
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed.vtm
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed_0_0.vtu
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed_1_0.vtu
      ${ALE_DATA_DIR}/meshes/2d/mixed/mixed.vtk
      ${ALE_DATA_DIR}/meshes/2d/voro/voro.g
      ${ALE_DATA_DIR}/meshes/3d/box/box-hex.g
      ${ALE_DATA_DIR}/meshes/3d/box/box-nfaced.exo
      ${ALE_DATA_DIR}/meshes/3d/box/box-tet.g
      ${ALE_DATA_DIR}/meshes/3d/voro/sedov-0_05cm.exo
    LIBRARIES ${mesh_LIBS} ale
)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
