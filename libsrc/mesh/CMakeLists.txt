#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(mesh_HEADERS
  burton/burton.h

  burton/burton_io_exodus.h
  burton/burton_io_tecplot.h
  burton/burton_io_vtk.h
  burton/burton_io_vtm.h
  burton/burton_io_vtu.h

  burton/burton_config.h
  burton/burton_mesh.h
  burton/burton_mesh_topology.h
  burton/burton_types.h

  burton/burton_corner.h
  burton/burton_element.h
  burton/burton_vertex.h
  burton/burton_wedge.h

  burton/burton_triangle.h
  burton/burton_quadrilateral.h
  burton/burton_polygon.h

  burton/burton_tetrahedron.h
  burton/burton_hexahedron.h
  burton/burton_polyhedron.h

  factory.h
  mesh_utils.h
  voronoi.h  detail/voronoi_impl.h
  vtk_utils.h

  PARENT_SCOPE
)


set(mesh_SOURCES
  burton/burton_element.cc
  burton/burton_hexahedron.cc
  burton/burton_polygon.cc
  burton/burton_polyhedron.cc
  burton/burton_quadrilateral.cc
  burton/burton_tetrahedron.cc
  burton/burton_triangle.cc
  burton/burton_wedge.cc
  burton/burton_vertex.cc
  PARENT_SCOPE
)


# get a relative path
file(RELATIVE_PATH DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ${ALE_DATA_DIR})

cinch_add_unit(mesh
    SOURCES 
      burton/test/burton_test_base.h
      burton/test/burton_create_test.h
      burton/test/burton_2d_test.h
      burton/test/burton_3d_test.h
      burton/test/burton_io_test.h
      burton/test/burton_voro_test.h

      burton/test/burton_create.cc
      burton/test/burton_2d.cc
      burton/test/burton_3d.cc
      burton/test/burton_io.cc
      burton/test/burton_voro.cc
      burton/test/burton_grad.cc
    INPUTS 
      burton/test/burton_2d.dump.blessed
      burton/test/burton_2d.mesh.blessed
      burton/test/burton_2d.geometry.blessed
      burton/test/burton_2d.copy.blessed
      burton/test/burton_3d.dump.blessed
      burton/test/burton_3d.mesh.blessed
      burton/test/burton_3d.geometry.blessed
      burton/test/burton_create.faces.blessed
      burton/test/burton_create.minimal.blessed
      burton/test/burton_create.points.blessed
      ${DATA_DIR}/meshes/2d/mixed/mixed.g
      ${DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${DATA_DIR}/meshes/2d/mixed/mixed.vtu
      ${DATA_DIR}/meshes/2d/mixed/mixed.vtm
      ${DATA_DIR}/meshes/2d/mixed/mixed_0_0.vtu
      ${DATA_DIR}/meshes/2d/mixed/mixed_1_0.vtu
      ${DATA_DIR}/meshes/2d/mixed/mixed.vtk
      ${DATA_DIR}/meshes/2d/voro/voro.g
      ${DATA_DIR}/meshes/3d/box/box-hex.exo
      ${DATA_DIR}/meshes/3d/box/box-hex.g
      ${DATA_DIR}/meshes/3d/box/box-nfaced.exo
      ${DATA_DIR}/meshes/3d/box/box-tet.exo
      ${DATA_DIR}/meshes/3d/box/box-tet.g
      ${DATA_DIR}/meshes/3d/voro/sedov-0_05cm.exo
    LIBRARIES ${ALE_LIBRARIES}
)
