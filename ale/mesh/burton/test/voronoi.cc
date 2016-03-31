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

//! system includes
#include<string>

//! user includes
#include "burton_test.h"
#include "../../../mesh/voronoi.h"

// some general using statements
using std::string;
using ale::mesh::write_mesh;
using ale::mesh::read_mesh;


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonVoro, constrained) {
  mesh_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));

  auto merge_tol = 1.e-12;
  auto num_steps = 10;
  auto use_clipping = false;
  auto voro = voronoi( m, merge_tol, num_steps, use_clipping );

  // write m to a different file
  name = "voro.vtk";
  ASSERT_FALSE(write_mesh(name, voro));

  // write m to a different file
  name = "voro.plt";
  ASSERT_FALSE(write_mesh(name, voro));

  // write m to a different file
  name = "voro.dat";
  ASSERT_FALSE(write_mesh(name, voro));

  // write m to a different file
  name = "voro.g";
  ASSERT_FALSE(write_mesh(name, voro));


} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST(BurtonVoro, read_write_exo) {
  mesh_t m;
  // read mesh written by above test
  string name("voro.g");
  ASSERT_FALSE(read_mesh(name, m));
  // write mesh written by above test
  name = "voro_out.g";
  ASSERT_FALSE(write_mesh(name, m));
} // TEST_F

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
