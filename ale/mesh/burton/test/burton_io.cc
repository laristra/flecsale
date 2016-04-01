/*~-------------------------------------------------------------------------~~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
/*!
 *
 * \file burton_io.cc
 * 
 * \brief Tests io of the burton mesh.
 *
 ******************************************************************************/

//! user includes
#include "burton_test.h"

// some general using statements
using std::string;
using ale::mesh::write_mesh;
using ale::mesh::read_mesh;

#ifdef HAVE_EXODUS

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, write_exo) {
  // create state data on b
  create_data();

  // write the mesh
  string name("mesh.exo");
  ASSERT_FALSE(write_mesh(name, mesh_));

  name = "mesh.g";
  ASSERT_FALSE(write_mesh(name, mesh_));

} // TEST_F

#endif // HAVE_EXODUS

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, write_dat) {
  // create state data on b
  create_data();
  // write the mesh
  string name("mesh.dat");
  ASSERT_FALSE(write_mesh(name, mesh_));
} // TEST_F



#ifdef HAVE_TECIO

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, write_plt) {
  // create state data on b
  create_data();
  // write the mesh
  string name("mesh.plt");
  ASSERT_FALSE(write_mesh(name, mesh_));
} // TEST_F

#endif // HAVE_TECIO

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, write_vtk) {
  // create state data on b
  create_data();
  // write the mesh in default format
  string name("mesh_default.vtk");
  ASSERT_FALSE(write_mesh(name, mesh_));
  // write the mesh in ascii
  name = "mesh_ascii.vtk";
  ASSERT_FALSE(write_mesh(name, mesh_, false));
  // write it again in binary
  name = "mesh_binary.vtk";
  ASSERT_FALSE(write_mesh(name, mesh_, true));
} // TEST_F


#ifdef HAVE_EXODUS

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonIO, read_exo) {
  mesh_t m;
  // read mesh written by above test
  string name("mesh.exo");
  ASSERT_FALSE(read_mesh(name, m));

  // write m to a different file
  name = "mesh_out.exo";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonIO, read_g) {
  mesh_t m;
  // read mesh written by above test
  string name("mesh.g");
  ASSERT_FALSE(read_mesh(name, m));

  // write m to a different file
  name = "mesh_out.g";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonIO, io_mixed) {
  mesh_t m;

  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));

  // write m to a different file
  name = "mixed_out.g";
  ASSERT_FALSE(write_mesh(name, m));

  // write m to a different file
  name = "mixed_out.vtk";
  ASSERT_FALSE(write_mesh(name, m));

#ifdef HAVE_VTK
  // write m to a different file
  name = "mixed_out.vtu";
  ASSERT_FALSE(write_mesh(name, m));

  // write m to a different file
  name = "mixed_out.vtm";
  ASSERT_FALSE(write_mesh(name, m));
#endif

#ifdef HAVE_TECIO
  // write m to a different file
  name = "mixed_out.plt";
  ASSERT_FALSE(write_mesh(name, m));
#endif

  // write m to a different file
  name = "mixed_out.dat";
  ASSERT_FALSE(write_mesh(name, m));

} // TEST_F


#endif // HAVE_EXODUS

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
