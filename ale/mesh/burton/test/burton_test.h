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
 * \file burton_test.h
 * 
 * \brief Defines a test fixture.
 *
 ******************************************************************************/
#pragma once

//! test include
#include "burton_test_base.h"

//! \brief the mesh float type
using real_t   = typename mesh_t::real_t;
//! \brief the mesh int type
using integer_t= typename mesh_t::integer_t;
//! \brief the mesh dimensions
static constexpr size_t dimensions = mesh_t::num_dimensions();

//! \brief the point
using point_t  = typename mesh_t::point_t;
//! \brief the vector type
using vector_t = typename mesh_t::vector_t;
//! \brief the vertex type
using vertex_t = typename mesh_t::vertex_t;
//! \brief the vertex type
using edge_t   = typename mesh_t::edge_t;
//! \brief the cell type
using cell_t   = typename mesh_t::cell_t;


// some general using statements
using std::size_t;
using std::vector;


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh
////////////////////////////////////////////////////////////////////////////////
class Burton : public BurtonTestBase {
public: 

  //---------------------------------------------------------------------------
  // Types
  //---------------------------------------------------------------------------
  
  //! \brief number of cells wide
  static constexpr size_t width = 2;
  //! \brief number of cells high
  static constexpr size_t height = 2;

  //! \brief some test tolerance
  static constexpr real_t test_tolerance = ale::common::test_tolerance;

protected:
  
  //---------------------------------------------------------------------------
  //! \brief the test setup function
  //! \remark this function is called before each test
  //---------------------------------------------------------------------------
  virtual void SetUp() {

    vector<vertex_t*> vs;
  
    mesh_.init_parameters((height+1)*(width+1));

    for(size_t j = 0; j < height + 1; ++j){
      for(size_t i = 0; i < width + 1; ++i){
        auto v = mesh_.create_vertex({ i, j });
        vs.push_back(v);
      } // for
    } // for

    auto  width1 = width + 1;

    for(size_t j = 0; j < height; ++j){
      for(size_t i = 0; i < width; ++i){
        // go over vertices counter clockwise to define cell
        auto c = mesh_.create_cell( {
            vs[i + j * width1],
            vs[i + 1 + j * width1],
            vs[i + 1 + (j + 1) * width1],
            vs[i + (j + 1) * width1] } );
      } // for
    } // for

    mesh_.init();

  } // SetUp

  //---------------------------------------------------------------------------
  //! \brief the test teardown function
  //! \remark this function is called after each test
  //---------------------------------------------------------------------------
  virtual void TearDown() { }



public:


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

  //! \brief the actual mesh object
  mesh_t mesh_;

};
