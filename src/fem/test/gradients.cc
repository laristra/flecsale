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

// system includes
#include <cinchtest.h>
#include <iostream>
#include <flexi/specializations/burton/burton.h>
#include <vector>


// explicitly use some stuff
using std::cout;
using std::endl;
using std::vector;
using flexi::persistent;

//=============================================================================
//! \brief Fixture for testing the partitioner.
//=============================================================================
class gradients : public ::testing::Test {

protected:

  //---------------------------------------------------------------------------
  // Types
  //---------------------------------------------------------------------------
  
  //! \brief the mesh type
  using mesh_t = flexi::burton_mesh_t;

  //! \brief the real type
  using real_t = mesh_t::real_t;

  //! \brief a dimensional space vector
  using vector_t = mesh_t::vector_t;

  //! \brief the vertex type
  using vertex_t = mesh_t::vertex_t;

  //! \brief the cell type
  using cell_t = mesh_t::cell_t;


  //---------------------------------------------------------------------------
  //! \brief the test setup function
  //---------------------------------------------------------------------------
  virtual void SetUp() {

    // reserve storage for the mesh
    mesh_.init_parameters((height+1)*(width+1));


    // create the individual vertices
    vector<vertex_t*> vs;
    
    for(size_t j = 0; j < height + 1; ++j){
      for(size_t i = 0; i < width + 1; ++i){
	auto v =
	  mesh_.create_vertex({double(i)+ 0.1*pow(double(j),1.8), 1.5*double(j)});
	v->set_rank(1);
	vs.push_back(v);
      }

    }

    // define each cell
    size_t width1 = width + 1;

    for(size_t j = 0; j < height; ++j){
      for(size_t i = 0; i < width; ++i){
	auto c = 
	  mesh_.create_cell({vs[i + j * width1],
		vs[i + 1 + j * width1],
		vs[i + 1 + (j + 1) * width1],
		vs[i + (j + 1) * width1]});
      }
    }

    // now finalize the mesh setup
    mesh_.init();
  }

  //---------------------------------------------------------------------------
  //! \brief the test teardown function 
  //---------------------------------------------------------------------------
  virtual void TearDown() { }


  //---------------------------------------------------------------------------
  //! \brief Determine the cell-cell adjacency information
  //!  
  //! \param[out] cell_idx   The cell-cell adjacency starting index for each cell
  //! \param[out] cell_neigh The cell-cell adjacency information
  //---------------------------------------------------------------------------
  
  template < typename E >
  decltype(auto) get_cell_neighbors( E cell ) {

    vector<E> neighbors;

    for ( auto e : mesh_.edges(cell) ) 
      for ( auto neigh : mesh_.cells(e) ) 
        if ( neigh != cell ) 
          neighbors.push_back( neigh );
      
    return neighbors;
  }

  
  //---------------------------------------------------------------------------
  // Data members
  //---------------------------------------------------------------------------

  //! \brief the mesh object used for testing
  mesh_t mesh_;

  //! \brief number of vertices in width of domain
  const size_t width = 10;
  //! \brief number of vertices in height of domain
  const size_t height = 20;

  //! \brief some test tolerance
  const real_t tol = 1e-12;

};


//=============================================================================
//! \brief Computation of simple cell centered gradients
//!
//! This test creates a simple flexi mesh, then computes gradients
//=============================================================================
TEST_F(gradients, simple_cell_centered) {

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
    auto x = mesh_.centroid(c);
    velocity[c] = x;
    std::cout << velocity[c][0] << " " << velocity[c][1] << std::endl;
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
      auto points = mesh_.vertices(e).toVec();
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
    auto x0 = mesh_.centroid(c0);
    auto u0 = const_velocity[c0];

    // get the neighbors
    auto neighbors = get_cell_neighbors(c0);

    // initialize the temporary variables
    real_t dxdx = 0;
    real_t dxdy = 0;
    real_t dydy = 0;
    decltype(u0) dudx = { 0, 0 };
    decltype(u0) dudy = { 0, 0 };

    // loop over each edge connected to v0
    for ( auto c1 : neighbors ) {

      // get the state and coordinates
      auto u1 = const_velocity[c1];
      auto x1 = mesh_.centroid(c1);

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
