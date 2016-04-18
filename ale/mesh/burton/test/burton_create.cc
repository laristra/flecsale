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
 * \file burton_create.cc
 * 
 * \brief Defines some mesh creation tests.
 *
 ******************************************************************************/

//! test include
#include "burton_create_test.h"


// some general using statements
using std::vector;


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh using the minimal set,
//!        i.e. a cell->vertex map only
////////////////////////////////////////////////////////////////////////////////
TEST_F( burton_create, minimal ) {

  // create mesh
  mesh_t mesh;

  // initialize  mesh
  mesh.init_parameters( num_vert );

  // create vertices
  vector<vertex_t*> vs( num_vert );

  for( size_t k = 0; k < num_vert_z; ++k )
    for( size_t j = 0; j < num_vert_y; ++j )
      for( size_t i = 0; i < num_vert_x; ++i )
      {
        auto v = mesh.create_vertex({ i, j, k });
        vs[ vert_index(i,j,k) ] =  v;
      }
  // create the cell
  for( size_t k = 0; k < num_cells_z; ++k )
    for( size_t j = 0; j < num_cells_y; ++j )
      for( size_t i = 0; i < num_cells_x; ++i )
        auto c = mesh.create_cell( 
          {
            vs[ vert_index( i  , j  , k  ) ],
            vs[ vert_index( i+1, j  , k  ) ],
            vs[ vert_index( i+1, j+1, k  ) ],
            vs[ vert_index( i  , j+1, k  ) ],
            vs[ vert_index( i  , j  , k+1) ],
            vs[ vert_index( i+1, j  , k+1) ],
            vs[ vert_index( i+1, j+1, k+1) ],
            vs[ vert_index( i  , j+1, k+1) ],
          } );

  // initialize the mesh
  mesh.init();

  // write m to a different file
  auto name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, mesh));

}

////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh using cell/face -> vertex
//!        maps only
////////////////////////////////////////////////////////////////////////////////
TEST_F( burton_create, points ) {

  // create mesh
  mesh_t mesh;

  // initialize  mesh
  mesh.init_parameters( num_vert );

  // create vertices
  vector<vertex_t*> vs( num_vert );

  for( size_t k = 0; k < num_vert_z; ++k )
    for( size_t j = 0; j < num_vert_y; ++j )
      for( size_t i = 0; i < num_vert_x; ++i )
      {
        auto v = mesh.create_vertex({ i, j, k });
        vs[ vert_index(i,j,k) ] =  v;
      }

  // create faces
  for( size_t k = 0; k < num_vert_z; ++k ) 
    for( size_t j = 0; j < num_vert_y; ++j ) 
      for( size_t i = 0; i < num_vert_x; ++i ) {

        // x plane
        if ( j < num_cells_y && k < num_cells_z )
          mesh.create_face( 
            {
              vs[ vert_index( i  , j  , k  ) ],
              vs[ vert_index( i  , j  , k+1) ],
              vs[ vert_index( i  , j+1, k+1) ],
              vs[ vert_index( i  , j+1, k  ) ]
            } );

        // y plane
        if ( i < num_cells_x && k < num_cells_z )
          mesh.create_face( 
            {
              vs[ vert_index( i  , j  , k  ) ],
              vs[ vert_index( i+1, j  , k  ) ],
              vs[ vert_index( i+1, j  , k+1) ],
              vs[ vert_index( i  , j  , k+1) ]
            } );

        // z plane
        if ( i < num_cells_x && j < num_cells_y )
          mesh.create_face( 
            {
              vs[ vert_index( i  , j  , k  ) ],
              vs[ vert_index( i+1, j  , k  ) ],
              vs[ vert_index( i+1, j+1, k  ) ],
              vs[ vert_index( i  , j+1, k  ) ]
            } );

      }
    
  // create the cell
  for( size_t k = 0; k < num_cells_z; ++k )
    for( size_t j = 0; j < num_cells_y; ++j )
      for( size_t i = 0; i < num_cells_x; ++i )
        auto c = mesh.create_cell( 
          {
            vs[ vert_index( i  , j  , k  ) ],
            vs[ vert_index( i+1, j  , k  ) ],
            vs[ vert_index( i+1, j+1, k  ) ],
            vs[ vert_index( i  , j+1, k  ) ],
            vs[ vert_index( i  , j  , k+1) ],
            vs[ vert_index( i+1, j  , k+1) ],
            vs[ vert_index( i+1, j+1, k+1) ],
            vs[ vert_index( i  , j+1, k+1) ],
          } );

  // initialize the mesh
  mesh.init();

  // write m to a different file
  auto name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, mesh));

}


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh using face->vertex and 
//!        cell->face maps
////////////////////////////////////////////////////////////////////////////////
TEST_F( burton_create, faces ) {

  // create mesh
  mesh_t mesh;

  // initialize  mesh
  mesh.init_parameters( num_vert );

  // create vertices
  vector<vertex_t*> vs( num_vert );

  for( size_t k = 0; k < num_vert_z; ++k )
    for( size_t j = 0; j < num_vert_y; ++j )
      for( size_t i = 0; i < num_vert_x; ++i )
      {
        auto v = mesh.create_vertex({ i, j, k });
        vs[ vert_index(i,j,k) ] =  v;
      }

  // storage for faces
  vector<face_t*> fs_x( num_face_x, nullptr );
  vector<face_t*> fs_y( num_face_y, nullptr );
  vector<face_t*> fs_z( num_face_z, nullptr );

  // create x-direction faces
  for( size_t k = 0; k < num_cells_z; ++k )
    for( size_t j = 0; j < num_cells_y; ++j ) 
      for( size_t i = 0; i < num_vert_x; ++i ) {

        auto f = mesh.create_face( 
            {
              vs[ vert_index( i  , j  , k  ) ],
              vs[ vert_index( i  , j  , k+1) ],
              vs[ vert_index( i  , j+1, k+1) ],
              vs[ vert_index( i  , j+1, k  ) ]
            } );         
        
        fs_x[ face_x_index(i,j,k) ] = f;
      }

  // create y-direction faces
  for( size_t k = 0; k < num_cells_z; ++k )
    for( size_t j = 0; j < num_vert_y; ++j ) 
      for( size_t i = 0; i < num_cells_x; ++i ) {

        auto f = mesh.create_face( 
            {
              vs[ vert_index( i  , j  , k  ) ],
              vs[ vert_index( i+1, j  , k  ) ],
              vs[ vert_index( i+1, j  , k+1) ],
              vs[ vert_index( i  , j  , k+1) ]
            } );         
        
        fs_y[ face_y_index(i,j,k) ] = f;
      }

  // create z-direction faces
  for( size_t k = 0; k < num_vert_z; ++k )
    for( size_t j = 0; j < num_cells_y; ++j ) 
      for( size_t i = 0; i < num_cells_x; ++i ) {

        auto f = mesh.create_face( 
            {
              vs[ vert_index( i  , j  , k  ) ],
              vs[ vert_index( i+1, j  , k  ) ],
              vs[ vert_index( i+1, j+1, k  ) ],
              vs[ vert_index( i  , j+1, k  ) ]
            } );         
        
        fs_z[ face_z_index(i,j,k) ] = f;
      }

  // create the cells
  for( size_t k = 0; k < num_cells_z; ++k )
    for( size_t j = 0; j < num_cells_y; ++j )
      for( size_t i = 0; i < num_cells_x; ++i ) 

        auto c = mesh.create_cell( 
          {
            fs_x[ face_x_index( i  , j  , k   ) ],
            fs_y[ face_y_index( i  , j  , k   ) ],
            fs_x[ face_x_index( i+1, j  , k   ) ],
            fs_y[ face_y_index( i  , j+1, k   ) ],
            fs_z[ face_z_index( i  , j  , k   ) ],
            fs_z[ face_z_index( i  , j  , k+1 ) ]
          } );

  // initialize the mesh
  mesh.init();

  // write m to a different file
  auto name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, mesh));

}

