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
 * \file gradients.cc
 * 
 * \brief Try computing some gradients on the burton mesh.
 *
 ******************************************************************************/

//! user includes
#include "burton_test.h"
#include "ale/geom/centroid.h"


// explicitly use some stuff
using std::cout;
using std::endl;
using ale::geom::centroid;

//=============================================================================
//! \brief Computation of simple cell centered gradients
//!
//! This test creates a simple flexi mesh, then computes gradients
//=============================================================================
TEST_F(Burton, gradients) {

  // get the number of cells in the mesh
  auto num_cells = mesh_.num_cells();
  auto num_verts = mesh_.num_vertices();

  //---------------------------------------------------------------------------
  // Create fields and assign gradients
  //---------------------------------------------------------------------------


  // setup a cell centered and vertex-based field
  register_state(mesh_, "pressure", vertices,   real_t, persistent);
  register_state(mesh_, "velocity",    cells, vector_t, persistent);

  // access the state arrays
  auto pressure = access_state(mesh_, "pressure",   real_t);
  auto velocity = access_state(mesh_, "velocity", vector_t);

  // loop over vertices and assign the field value
  for ( auto v : mesh_.vertices() ) {
    auto x = v->coordinates();
    pressure[v] = x[0] + x[1];
  }

  // loop over vertices and assign the field value
  for ( auto c : mesh_.cells() ) {
    auto x = centroid(c);
    velocity[c] = x;
  }

  //---------------------------------------------------------------------------
  // Compute gradients of point-based data
  //---------------------------------------------------------------------------


  // get constant accessors to the pressure
  const auto const_pressure = access_state(mesh_, "pressure",   real_t);

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
      auto points = mesh_.vertices(e).to_vec();
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
    ASSERT_NEAR( 1, dudx_ave, tol );
    ASSERT_NEAR( 1, dudy_ave, tol );
  }

  //---------------------------------------------------------------------------
  // Compute gradients of point-based data
  //---------------------------------------------------------------------------


  // get constant accessors to the pressure
  const auto const_velocity = access_state(mesh_, "velocity", vector_t);

  // loop over vertices and compute the gradient
  for ( auto c0 : mesh_.cells() ) {

    // the the current point and state
    auto x0 = centroid(c0);
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
      auto x1 = centroid(c1);

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
    ASSERT_NEAR( 1, dudx_ave[0], tol );
    ASSERT_NEAR( 0, dudx_ave[1], tol );

    ASSERT_NEAR( 0, dudy_ave[0], tol );
    ASSERT_NEAR( 1, dudy_ave[1], tol );
  }



  //---------------------------------------------------------------------------
  // Final Checks
  //---------------------------------------------------------------------------

  //ASSERT_TRUE(CINCH_EQUAL_BLESSED("scotch.blessed"));
  
    
} // TEST_F


/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
