/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Defines a test fixture.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// test include
#include "burton_test_base.h"

// user includes
#include "ale/mesh/factory.h"

// some general using statements
using std::vector;


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh
////////////////////////////////////////////////////////////////////////////////
class burton_2d : public burton_test_base {
public: 

  //---------------------------------------------------------------------------
  // Types From Mesh
  //---------------------------------------------------------------------------

  //! \brief the mesh type
  using mesh_t = mesh_2d_t;

  //! \brief the size type
  using size_t= typename mesh_t::size_t;
  //! \brief the counter type
  using counter_t= typename mesh_t::counter_t;
  //! \brief the mesh float type
  using real_t   = typename mesh_t::real_t;
  //! \brief the mesh int type
  using integer_t= typename mesh_t::integer_t;
  //! \brief the mesh dimensions
  static constexpr size_t num_dimensions = mesh_t::num_dimensions;

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

    mesh_ = ale::mesh::box<mesh_t>( 
      width, height, 0, 0, width, height );

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
