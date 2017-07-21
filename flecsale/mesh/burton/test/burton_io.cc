/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Tests io of the burton mesh.
///
////////////////////////////////////////////////////////////////////////////////

// user includes
#include "burton_io_test.h"

// register data
flecsi_register_data(mesh, hydro, pressure, real_t, dense, 1, cells);
flecsi_register_data(mesh, hydro, region, integer_t, dense, 1, cells);
flecsi_register_data(mesh, hydro, velocity, vector_t, dense, 1, vertices);


// Below tests need exodus to read the file
#ifdef HAVE_EXODUS 

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_exo_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.exo");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_exo_3d_poly) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-nfaced.exo");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_exo_3d_hex) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-hex.exo");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_exo_3d_tet) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-tet.exo");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_g_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".g";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_g_3d) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-hex.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".g";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_voro_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("voro.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // write mesh written by above test
  name = output_prefix()+".g";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, write_dat_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh
  name = output_prefix()+".dat";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F



////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an tecplot file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, write_dat_3d) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-tet.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".dat";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

#ifdef HAVE_TECIO

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, write_plt_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh
  name = output_prefix()+".plt";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading/writing  an tecplot file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, write_plt_3d) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-tet.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write m to a different file
  name = output_prefix()+".plt";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F


#endif // HAVE_TECIO

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing a vtk file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, write_vtk_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh in default format
  name = output_prefix()+"-default.vtk";
  ASSERT_FALSE(write_mesh(name, m));
  // write the mesh in ascii
  name = output_prefix()+"-ascii.vtk";
  ASSERT_FALSE(write_mesh(name, m, false));
  // write it again in binary
  name = output_prefix()+"-binary.vtk";
  ASSERT_FALSE(write_mesh(name, m, true));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing a vtk file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, write_vtk_3d) {
  mesh_3d_t m;
  // read mesh written by above test
  string name("box-tet.g");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh in default format
  name = output_prefix()+"-default.vtk";
  ASSERT_FALSE(write_mesh(name, m));
  // write the mesh in ascii
  name = output_prefix()+"-ascii.vtk";
  ASSERT_FALSE(write_mesh(name, m, false));
  // write it again in binary
  name = output_prefix()+"-binary.vtk";
  ASSERT_FALSE(write_mesh(name, m, true));
} // TEST_F

#endif // HAVE_EXODUS


// Below tests have their own readers

#ifdef HAVE_VTK

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an vtk file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_vtk_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.vtk");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh in default format
  name = output_prefix()+".vtk";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_vtu_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.vtu");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh in default format
  name = output_prefix()+".vtu";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_io, read_write_vtm_2d) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.vtm");
  ASSERT_FALSE(read_mesh(name, m));
  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );
  // create state data on b
  create_data(m);
  // write the mesh in default format
  name = output_prefix()+".vtm";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

#endif // HAVE_VTK
