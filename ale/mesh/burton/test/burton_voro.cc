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

//! test includes
#include "burton_voro_test.h"

//! user includes
#include "../../../mesh/voronoi.h"


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(BurtonVoro, constrained) {
  mesh_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));

  auto merge_tol = 1.e-12;
  auto num_steps = 10;
  auto use_clipping = false;
  auto voro = voronoi( m, merge_tol, num_steps, use_clipping );

  // write m to a different file
  name = output_prefix()+".vtk";
  ASSERT_FALSE(write_mesh(name, voro));

  // write m to a different file
  name = output_prefix()+".plt";
  ASSERT_FALSE(write_mesh(name, voro));

  // write m to a different file
  name = output_prefix()+".dat";
  ASSERT_FALSE(write_mesh(name, voro));

  // write m to a different file
  name = output_prefix()+".g";
  ASSERT_FALSE(write_mesh(name, voro));


} // TEST_F


/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
