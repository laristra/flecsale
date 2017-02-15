/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Try computing some gradients on the burton mesh.
///
////////////////////////////////////////////////////////////////////////////////

// user includes
#include "burton_2d_test.h"


// explicitly use some stuff
using std::cout;
using std::endl;

//=============================================================================
//! \brief Computation of simple cell centered gradients
//!
//! This test creates a simple flexi mesh, then computes gradients
//=============================================================================
TEST_F(burton_2d, gradients) {

  // get the number of cells in the mesh
  auto num_cells = mesh_.num_cells();
  auto num_verts = mesh_.num_vertices();

  //---------------------------------------------------------------------------
  // Create fields and assign gradients
  //---------------------------------------------------------------------------


  // setup a cell centered and vertex-based field
  register_data(mesh_, hydro, pressure, real_t, dense, 1, vertices);
  register_data(mesh_, hydro, velocity, vector_t, dense, 1, cells);

  // access the state arrays
  auto pressure = get_accessor(mesh_, hydro, pressure,   real_t, dense, 0);
  auto velocity = get_accessor(mesh_, hydro, velocity, vector_t, dense, 0);

  // loop over vertices and assign the field value
  for ( auto v : mesh_.vertices() ) {
    auto x = v->coordinates();
    pressure[v] = x[0] + x[1];
  }

  // loop over vertices and assign the field value
  for ( auto c : mesh_.cells() ) {
    auto x = c->centroid();
    velocity[c] = x;
  }

  //---------------------------------------------------------------------------
  // Compute gradients of point-based data
  //---------------------------------------------------------------------------


  // get constant accessors to the pressure
  const auto const_pressure = 
    get_accessor(mesh_, hydro, pressure, real_t, dense, 0);

  // loop over vertices and compute the gradient
  for ( auto v0 : mesh_.vertices() ) {

    // the the current point and state
    auto x0 = v0->coordinates();
    auto u0 = const_pressure[v0];

    // initialize the temporary variables
    real_t dxdx = 0;
    real_t dxdy = 0;
    real_t dydy = 0;
    decltype(u0) dudx = 0;
    decltype(u0) dudy = 0;

    // loop over each edge connected to v0
    for ( auto e : mesh_.edges(v0) ) {
      // get the edge points and select the other point
      auto points = mesh_.vertices(e);
      auto v1 = ( points[0] == v0 ? points[1] : points[0] );
      // get the state and coordinates
      auto u1 = const_pressure[v1];
      auto x1 = v1->coordinates();
      // compute the terms necessary for gradients
      auto dx = x1[0]-x0[0];
      auto dy = x1[1]-x0[1];
      dxdx += dx*dx;
      dxdy += dx*dy;
      dydy += dy*dy;
      auto du = u1-u0;
      dudx += du*dx;
      dudy += du*dy;
    }
    
    // the final gradient computation
    auto denom = 1 / ( dxdx*dydy-dxdy*dxdy );
    auto dudx_ave = 
      (dudx*dydy-dudy*dxdy) * denom;
    auto dudy_ave = 
      (dudy*dxdx-dudx*dxdy) * denom;

    // make sure the gradient is 1
    ASSERT_NEAR( 1, dudx_ave, test_tolerance );
    ASSERT_NEAR( 1, dudy_ave, test_tolerance );
  }

  //---------------------------------------------------------------------------
  // Compute gradients of point-based data
  //---------------------------------------------------------------------------


  // get constant accessors to the pressure
  const auto const_velocity = 
    get_accessor(mesh_, hydro, velocity, vector_t, dense, 0);

  // loop over vertices and compute the gradient
  for ( auto c0 : mesh_.cells() ) {

    // the the current point and state
    auto x0 = c0->centroid();
    auto u0 = const_velocity[c0];

    // get the neighbors
    auto neighbors = get_cell_neighbors(c0);

    // initialize the temporary variables
    real_t dxdx = 0;
    real_t dxdy = 0;
    real_t dydy = 0;
    decltype(u0) dudx = { 0., 0. };
    decltype(u0) dudy = { 0., 0. };

    // loop over each edge connected to v0
    for ( auto c1 : neighbors ) {

      // get the state and coordinates
      auto u1 = const_velocity[c1];
      auto x1 = c1->centroid();

      // compute the terms necessary for gradients
      auto dx = x1[0]-x0[0];
      auto dy = x1[1]-x0[1];
      dxdx += dx*dx;
      dxdy += dx*dy;
      dydy += dy*dy;
      auto du = u1-u0;
      dudx += du*dx;
      dudy += du*dy;
    }
    
    // the final gradient computation
    auto denom = 1 / ( dxdx*dydy-dxdy*dxdy );
    auto dudx_ave = 
      (dudx*dydy-dudy*dxdy) * denom;
    auto dudy_ave = 
      (dudy*dxdx-dudx*dxdy) * denom;

    // make sure the gradient is du/dx = {1,0} and du/dy = {0,1}
    ASSERT_NEAR( 1, dudx_ave[0], test_tolerance );
    ASSERT_NEAR( 0, dudx_ave[1], test_tolerance );

    ASSERT_NEAR( 0, dudy_ave[0], test_tolerance );
    ASSERT_NEAR( 1, dudy_ave[1], test_tolerance );
  }



  //---------------------------------------------------------------------------
  // Final Checks
  //---------------------------------------------------------------------------

  //ASSERT_TRUE(CINCH_EQUAL_BLESSED("scotch.blessed"));
  
    
} // TEST_F

