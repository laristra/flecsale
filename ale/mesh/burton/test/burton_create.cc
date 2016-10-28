/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Defines some mesh creation tests.
///
////////////////////////////////////////////////////////////////////////////////

// test include
#include "burton_create_test.h"


// some general using statements
using std::endl;
using std::vector;

////////////////////////////////////////////////////////////////////////////////
//! \brief A utility to dump the connectivity of a 3d mesh.
////////////////////////////////////////////////////////////////////////////////
template< typename T >
auto & dump_connectivity( std::ostream & os, T & mesh ) 
{
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  os << separator;
  os << "mesh" << endl;

  os << separator;
  os << "Vertices in mesh:" << endl;

  for(auto v : mesh.vertices())
    os << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;

  os << separator;
  os << "Edges in mesh:" << endl;

  for(auto e : mesh.edges())
    os << "----------- edge id: " << e.id()
                    << " with midpoint " << e->midpoint() << endl;

  os << separator;
  os << "Faces in mesh:" << endl;

  for(auto f : mesh.faces()) 
    os << "----------- faces id: " << f.id()
                    << " with centroid " << f->centroid() << endl;

  os << separator;
  os << "Cells in mesh:" << endl;

  for(auto c : mesh.cells()) 
    os << "----------- cell id: " << c.id()
                    << " with centroid " << c->centroid() << endl;

  os << separator;
  os << "Corners in mesh:" << endl;

  for(auto c : mesh.corners()) {
    os << "----------- corner id: " << c.id() << endl;
  } // for

  os << separator;
  os << "Wedges in mesh:" << endl;

  for(auto w : mesh.wedges()) {
    os << "----------- wedge id: " << w.id() << endl;
  } // for

  os << separator;
  os << "For each vertex:" << endl;

  for(auto v: mesh.vertices()) {
    os << "^^^^^^^^Vertex id: " << v.id() << endl;

    os << "    ----Corners:" << endl;
    for(auto c: mesh.corners(v))
      os << "    ++++ corner id: " << c.id() << endl;

    os << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(v))
      os << "    ++++ wedge id: " << w.id() << endl;

    os << "    ----Cells:" << endl;
    for(auto c: mesh.cells(v))
      os << "    ++++ cell id: " << c.id() << endl;

    os << "    ----Faces:" << endl;
    for(auto f: mesh.faces(v))
      os << "    ++++ face id: " << f.id() << endl;

    os << "    ----Edges:" << endl;
    for(auto e: mesh.edges(v))
      os << "    ++++ edge id: " << e.id() << endl;
  } // for

  os << separator;
  os << "For each edge:" << endl;

  for(auto e : mesh.edges()) {
    os << "^^^^^^^^Edge id: " << e.id() << endl;

    os << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(e))
      os << "    ++++ corner id: " << cnr.id() << endl;

    os << "    ----Cells:" << endl;
    for(auto c : mesh.cells(e))
      os << "    ++++ cell id: " << c.id() << endl;

    os << "    ----Faces:" << endl;
    for(auto f: mesh.faces(e))
      os << "    ++++ face id: " << f.id() << endl;

    os << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(e))
      os << "    ++++ vertex id: " << v.id() << endl;

  } // for

  os << separator;
  os << "For each face:" << endl;

  for(auto f : mesh.faces()) {
    os << "^^^^^^^^Face id: " << f.id() << endl;

    os << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(f))
      os << "    ++++ corner id: " << cnr.id() << endl;

    os << "    ----Cells:" << endl;
    for(auto c: mesh.cells(f))
      os << "    ++++ cell id: " << c.id() << endl;

    os << "    ----Edges:" << endl;
    for(auto e : mesh.edges(f))
      os << "    ++++ edge id: " << e.id() << endl;

    os << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(f))
      os << "    ++++ vertex id: " << v.id() << endl;

  } // for

  os << separator;
  os << "For each cell:" << endl;

  for(auto c : mesh.cells()) {
    os << "^^^^^^^^Cell id: " << c.id() << endl;

    os << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(c))
      os << "    ++++ corner id: " << cnr.id() << endl;

    os << "    ----Faces:" << endl;
    for(auto f: mesh.faces(c))
      os << "    ++++ face id: " << f.id() << endl;

    os << "    ----Edges:" << endl;
    for(auto e : mesh.edges(c))
      os << "    ++++ edge id: " << e.id() << endl;

    os << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(c))
      os << "    ++++ vertex id: " << v.id() << endl;

  } // for


  os << separator;
  os << "For each corner:" << endl;

  for(auto c : mesh.corners()) {
    os << "^^^^^^^^Corner id: " << c.id() << endl;

    os << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(c)) 
      os << "    ++++ wedge id: " << w.id() << endl;

    os << "    ----Cells:" << endl;
    for(auto cl: mesh.cells(c)) 
      os << "    ++++ cell id: " << cl.id() << endl;

    os << "    ----Faces:" << endl;
    for(auto f: mesh.faces(c))
      os << "    ++++ face id: " << f.id() << endl;

    os << "    ----Edges:" << endl;
    for(auto e: mesh.edges(c)) 
      os << "    ++++ edge id: " << e.id() << endl;

    os << "    ----Vertices:" << endl;
    for(auto v: mesh.vertices(c)) 
      os << "    ++++ vertex id: " << v.id() << endl;

  } // for

  os << separator;
  os << "For each wedge:" << endl;

  for(auto w : mesh.wedges()) {
    os << "^^^^^^^^Wedge id: " << w.id() << endl;

    os << "    ----Corners:" << endl;
    for(auto c: mesh.corners(w)) 
      os << "    ++++ corner id: " << c.id() << endl;

    os << "    ----Cells:" << endl;
    for(auto cl: mesh.cells(w)) 
      os << "    ++++ cell id: " << cl.id() << endl;

    os << "    ----Faces:" << endl;
    for(auto f: mesh.faces(w))
      os << "    ++++ face id: " << f.id() << endl;

    os << "    ----Edges:" << endl;
    for(auto e: mesh.edges(w)) 
      os << "    ++++ edge id: " << e.id() << endl;

    os << "    ----Vertices:" << endl;
    for(auto v: mesh.vertices(w)) 
      os << "    ++++ vertex id: " << v.id() << endl;

  } // for

  // return the stream
  return os;
}

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

  for( counter_t k = 0; k < num_vert_z; ++k )
    for( counter_t j = 0; j < num_vert_y; ++j )
      for( counter_t i = 0; i < num_vert_x; ++i )
      {
        auto v = mesh.create_vertex({ i, j, k });
        vs[ vert_index(i,j,k) ] =  v;
      }
  // create the cell
  for( counter_t k = 0; k < num_cells_z; ++k )
    for( counter_t j = 0; j < num_cells_y; ++j )
      for( counter_t i = 0; i < num_cells_x; ++i )
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

  // dump and check the connectivity
  dump_connectivity( CINCH_CAPTURE(), mesh );
  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

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

  for( counter_t k = 0; k < num_vert_z; ++k )
    for( counter_t j = 0; j < num_vert_y; ++j )
      for( counter_t i = 0; i < num_vert_x; ++i )
      {
        auto v = mesh.create_vertex({ i, j, k });
        vs[ vert_index(i,j,k) ] =  v;
      }

  // create faces
  for( counter_t k = 0; k < num_vert_z; ++k ) 
    for( counter_t j = 0; j < num_vert_y; ++j ) 
      for( counter_t i = 0; i < num_vert_x; ++i ) {

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
  for( counter_t k = 0; k < num_cells_z; ++k )
    for( counter_t j = 0; j < num_cells_y; ++j )
      for( counter_t i = 0; i < num_cells_x; ++i )
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

  // dump and check the connectivity
  dump_connectivity( CINCH_CAPTURE(), mesh );
  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

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

  for( counter_t k = 0; k < num_vert_z; ++k )
    for( counter_t j = 0; j < num_vert_y; ++j )
      for( counter_t i = 0; i < num_vert_x; ++i )
      {
        auto v = mesh.create_vertex({ i, j, k });
        vs[ vert_index(i,j,k) ] =  v;
      }

  // storage for faces
  vector<face_t*> fs_x( num_face_x, nullptr );
  vector<face_t*> fs_y( num_face_y, nullptr );
  vector<face_t*> fs_z( num_face_z, nullptr );

  // create x-direction faces
  for( counter_t k = 0; k < num_cells_z; ++k )
    for( counter_t j = 0; j < num_cells_y; ++j ) 
      for( counter_t i = 0; i < num_vert_x; ++i ) {

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
  for( counter_t k = 0; k < num_cells_z; ++k )
    for( counter_t j = 0; j < num_vert_y; ++j ) 
      for( counter_t i = 0; i < num_cells_x; ++i ) {

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
  for( counter_t k = 0; k < num_vert_z; ++k )
    for( counter_t j = 0; j < num_cells_y; ++j ) 
      for( counter_t i = 0; i < num_cells_x; ++i ) {

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
  for( counter_t k = 0; k < num_cells_z; ++k )
    for( counter_t j = 0; j < num_cells_y; ++j )
      for( counter_t i = 0; i < num_cells_x; ++i ) 

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

  // dump and check the connectivity
  dump_connectivity( CINCH_CAPTURE(), mesh );
  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

  // write m to a different file
  auto name = output_prefix()+".exo";
  ASSERT_FALSE(write_mesh(name, mesh));

}

