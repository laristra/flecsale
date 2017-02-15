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

//! test includes
#include "burton_voro_test.h"

// user includes
#include "mesh/voronoi.h"

#if defined(HAVE_EXODUS) && defined(HAVE_SHAPO)

////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_voro, constrained) {
  mesh_2d_t m;
  // read mesh written by above test
  string name("mixed.g");
  ASSERT_FALSE(read_mesh(name, m));

  auto merge_tol = 1.e-12;
  auto num_steps = 0;
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


////////////////////////////////////////////////////////////////////////////////
//! \brief test reading an exodus file
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_voro, read_write_3d) {
  mesh_3d_t m;

  // read mesh written by above test
  string name("sedov-0_05cm.exo");
  ASSERT_FALSE(read_mesh(name, m));

  // check the mesh
  EXPECT_TRUE( m.is_valid(false) );

  // write m to a different file
  name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, m));


} // TEST_F

#endif
