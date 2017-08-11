#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(mesh_HEADERS
  burton/burton.h

  burton/burton_config.h
  burton/burton_mesh.h
  burton/burton_mesh_topology.h
  burton/burton_types.h

  burton/burton_corner.h
  burton/burton_element.h
  burton/burton_vertex.h
  burton/burton_wedge.h

  factory.h
  mesh_utils.h

  portage/portage.h
  portage/portage_mesh.h
  portage/portage_state.h

  voronoi.h  detail/voronoi_impl.h
  vtk_utils.h

)


set(mesh_SOURCES
  burton/burton_element.cc
  burton/burton_types.cc
  burton/burton_wedge.cc
)

cinch_install_headers(
  FILES ${mesh_HEADERS}
  DESTINATION include/flecsale/mesh
)

add_library( flecsale_mesh OBJECT ${mesh_SOURCES} )
set(FleCSALE_OBJECTS ${FleCSALE_OBJECTS} flecsale_mesh PARENT_SCOPE)

# get a relative path
cinch_add_unit(test_mesh
    SOURCES 
      burton/test/burton_create.cc
      burton/test/burton_2d.cc
      burton/test/burton_3d.cc
      burton/test/burton_io.cc
      burton/test/burton_voro.cc
      burton/test/burton_grad.cc

      portage/test/portage_2d.cc
      portage/test/portage_3d.cc

      $<TARGET_OBJECTS:flecsale_mesh>
      $<TARGET_OBJECTS:flecsale_io>
      
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
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed.g
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed.vtu
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed.vtm
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed_0_0.vtu
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed_1_0.vtu
      ${FleCSALE_DATA_DIR}/meshes/2d/mixed/mixed.vtk
      ${FleCSALE_DATA_DIR}/meshes/2d/voro/voro.g
      ${FleCSALE_DATA_DIR}/meshes/3d/box/box-hex.exo
      ${FleCSALE_DATA_DIR}/meshes/3d/box/box-hex.g
      ${FleCSALE_DATA_DIR}/meshes/3d/box/box-nfaced.exo
      ${FleCSALE_DATA_DIR}/meshes/3d/box/box-tet.exo
      ${FleCSALE_DATA_DIR}/meshes/3d/box/box-tet.g
      ${FleCSALE_DATA_DIR}/meshes/3d/voro/sedov-0_05cm.exo
)

# make sure the target was actually made
if (TARGET test_mesh)
  # the libraries used for this sub-module.
  # link them all because there are quite a few of them.
  target_link_libraries(test_mesh ${FleCSALE_LIBRARIES})
endif()

