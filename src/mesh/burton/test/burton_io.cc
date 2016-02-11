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
using flecsi::write_mesh;
using flecsi::read_mesh;

////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, write_exo) {
  // create state data on b
  // register
  register_state(mesh_, "pressure", cells, real_t, persistent);
  register_state(mesh_, "region", cells, int, persistent);
  register_state(mesh_, "velocity", vertices, vector_t, persistent);
  // access
  auto p = access_state(mesh_, "pressure", real_t);
  auto r = access_state(mesh_, "region", int);
  auto velocity = access_state(mesh_, "velocity", vector_t);
  // initialize
  for(auto c: mesh_.cells()) {
    p[c] = c.id();
    r[c] = mesh_.num_cells() - c.id();
  } // for
  // vertices
  for (auto v: mesh_.vertices()) {
    velocity[v][0] = v.id();
    velocity[v][1] = 2.0*v.id();
  } // for

  // write the mesh
  string name("test/mesh.exo");
  ASSERT_FALSE(write_mesh(name, mesh_));

} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test writing an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, write_g) {
  string name("test/mesh.g");
  ASSERT_FALSE(write_mesh(name, mesh_));
} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading  an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonIO, read_exo) {
  mesh_t m;
  // read mesh written by above test
  string name("test/mesh.exo");
  ASSERT_FALSE(read_mesh(name, m));

  // write m to a different file
  name = "test/mesh_out.exo";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonIO, read_g) {
  mesh_t m;
  // read mesh written by above test
  string name("test/mesh.g");
  ASSERT_FALSE(read_mesh(name, m));

  // write m to a different file
  name = "test/mesh_out.g";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
