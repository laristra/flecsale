/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
/*!
 * \file burton_entity_types.h
 * \authors bergen
 * \date Initial file creation: Nov 15, 2015
 ******************************************************************************/

#pragma once

//! user includes
#include "ale/geom/shapes/polyhedron.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_polyhedron_t burton_entity_types.h
//!
//! \brief The burton_polyhedron_t type provides a derived instance of
//!   burton_cell_t for 2D polyhedron cells.
////////////////////////////////////////////////////////////////////////////////
class burton_polyhedron_t : public burton_element_t<3,3>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<3,3>;

  //============================================================================
  // Constructors
  //============================================================================

  //! main constructor
  template< typename F >
  burton_polyhedron_t(mesh_topology_base_t & mesh, F && faces) : base_t(mesh)
  {
    auto num_faces = std::forward<F>(faces).size();
    
  }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  point_t centroid() const override;

  //! the midpoint
  point_t midpoint() const override;

  //! the area of the cell
  real_t volume() const override;

  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::polyhedron<point_t>::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_polyhedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_entities(
    size_t dim, const id_t & cell,
    connectivity_t*  (&conn)[num_dimensions+1][num_dimensions], 
    id_t * entities )  override
  {        

    // you should only be coming in here to build edges
    assert( dim == 1 ); 

    // get the cell id
    auto cell_id = cell.entity();
    // get the cell entities
    size_t num_cell_faces = 0;
    auto cell_faces_begin = conn[3][2]->get_entities( cell_id, num_cell_faces );
    auto cell_faces_end   = std::next( cell_faces_begin, num_cell_faces );
    // make sure the faces exist
    assert( num_cell_faces > 0 && "no cell faces yet" );
    
    // the list of edge pairs
    std::vector< std::pair< id_t, id_t > > cell_edges;

    // get the edges of each face
    for ( auto face = cell_faces_begin;  face != cell_faces_end; face++ ) { 
      // get the face id
      auto face_id = face->entity();
      // get all vertices in this face
      size_t num_face_verts = 0;
      auto face_verts_begin = conn[2][0]->get_entities( face_id, num_face_verts ); 
      auto face_verts_end   = std::next( face_verts_begin, num_face_verts );      
      assert( num_face_verts > 2 && "not enough vertices for a valid face" );
      // reserve space 
      cell_edges.reserve( cell_edges.size() + num_face_verts );
      // add each edge pair to the list
      auto vp = std::prev( face_verts_end );
      assert( vp != face_verts_begin && "no vertices in this face" ); 
      for ( auto vn = face_verts_begin; vn != face_verts_end; vn++ ) {
        assert( *vn != *vp && "edge has two equal vertices" );
        if ( *vn < *vp ) 
          cell_edges.emplace_back( std::make_pair( *vn, *vp ) );
        else
          cell_edges.emplace_back( std::make_pair( *vp, *vn ) );
        vp = vn;
      }          
    }

    // now sort the list of edges
    std::sort( 
      cell_edges.begin(), cell_edges.end(), 
      [](const auto & a, const auto & b) 
      { 
        size_t a1 = a.first;
        size_t b1 = b.first;
        if ( a1 == b1 )
          return ( a.second < b.second );
        else
          return (a1 < b1);
      }
    );
    // remove uniques
    auto end = std::unique( cell_edges.begin(), cell_edges.end() );
    auto num_edges = std::distance( cell_edges.begin(), end );

    // copy the unique results to the output array   
    size_t i = 0;

    std::for_each( 
      cell_edges.begin(), end, 
      [&](const auto & edge) 
      {
        entities[i++] = edge.first;
        entities[i++] = edge.second;
      }
    );

    return std::vector<id_t>( num_edges, 2 );

      
  } // create_entities

  //----------------------------------------------------------------------------
  // \brief create_bound_entities function for burton_polyhedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    connectivity_t*  (&from_domain_conn)[num_dimensions+1][num_dimensions+1], 
    connectivity_t*  (&  to_domain_conn)[num_dimensions+1][num_dimensions+1], 
    id_t * entities )  override
  {
    // get the cell id
    auto cell_id = cell.entity();

    size_t i = 0;

    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      //
      // Take your right hand, its origin is the vertex of the corner.  Curl 
      // your hand from the first edge to the second edge, with the third edge
      // aligned with your thumb.  You hand also curls from the first to the 
      // first to second face, with the third face on the bottom.
      //
    case 0: {

      std::vector<id_t> entity_count;

      // get the cell entities
      size_t num_cell_verts = 0, num_cell_edges = 0, num_cell_faces = 0;
      auto cell_verts_begin = from_domain_conn[3][0]->get_entities( cell_id, num_cell_verts );
      auto cell_edges_begin = from_domain_conn[3][1]->get_entities( cell_id, num_cell_edges );
      auto cell_faces_begin = from_domain_conn[3][2]->get_entities( cell_id, num_cell_faces );
      auto cell_verts_end   = std::next( cell_verts_begin, num_cell_verts );
      auto cell_edges_end   = std::next( cell_edges_begin, num_cell_edges );
      auto cell_faces_end   = std::next( cell_faces_begin, num_cell_faces );

      // sort the edges and faces for intersections later      
      auto cell_edges = std::vector<id_t>( cell_edges_begin, cell_edges_end );
      auto cell_faces = std::vector<id_t>( cell_faces_begin, cell_faces_end );
      std::sort( cell_edges.begin(), cell_edges.end() );
      std::sort( cell_faces.begin(), cell_faces.end() );

      // temparary lists
      std::vector<id_t> edges, faces;

      for ( auto vert = cell_verts_begin;  vert != cell_verts_end; vert++ ) { 
        // get the vertex id
        auto vert_id = vert->entity();
        // clear temporary lits
        edges.clear();
        faces.clear();
        // get all entities attached to this vertex
        size_t num_vert_edges = 0, num_vert_faces = 0;
        auto vert_edges_begin = from_domain_conn[0][1]->get_entities( vert_id, num_vert_edges ); 
        auto vert_faces_begin = from_domain_conn[0][2]->get_entities( vert_id, num_vert_faces ); 
        auto vert_edges = std::vector<id_t>( vert_edges_begin, std::next( vert_edges_begin, num_vert_edges ) );
        auto vert_faces = std::vector<id_t>( vert_faces_begin, std::next( vert_faces_begin, num_vert_faces ) );
        // sort the lists for intersectinos
        std::sort( vert_edges.begin(), vert_edges.end() );
        std::sort( vert_faces.begin(), vert_faces.end() );
        // get the intersections of the sets
        std::set_intersection( vert_edges.begin(), vert_edges.end(),
                               cell_edges.begin(), cell_edges.end(),
                               std::back_inserter(edges));
        std::set_intersection( vert_faces.begin(), vert_faces.end(),
                               cell_faces.begin(), cell_faces.end(),
                               std::back_inserter(faces));
        // add the entities to the list
        entities[i++] = *vert;
        for ( auto e : edges ) entities[i++] = e;
        for ( auto f : faces ) entities[i++] = f;
        // add the final number of entities to the list
        entity_count.emplace_back( 1 + edges.size() + faces.size() );
        
      }  // foreach vertex

      return entity_count;
    }
      //------------------------------------------------------------------------
      // Wedges
      //
      // Each corner has 6 vertices.  There is an even/odd ordering so we 
      // no which way to compute normals.  So edges are defined in pairs
      //
    case 1: {
  
      // get the higher dimensional entities
      size_t num_cell_faces = 0, num_cell_corners = 0;
      auto cell_faces_begin = from_domain_conn[3][2]->get_entities( cell_id, num_cell_faces );
      auto cell_corners_begin = to_domain_conn[3][0]->get_entities( cell_id, num_cell_corners );
      auto cell_faces_end   = std::next( cell_faces_begin, num_cell_faces );
      auto cell_corners_end = std::next( cell_corners_begin, num_cell_corners );
      
      // a counter for the number of wedges
      size_t num_wedges = 0;

      // loop over faces
      for ( auto face = cell_faces_begin;  face != cell_faces_end; face++ ) { 

        // get the face id
        auto face_id = face->entity();

        // get the vertices of the face
        size_t num_face_verts = 0;
        auto face_verts_begin = from_domain_conn[2][0]->get_entities( face_id, num_face_verts );
        auto face_verts_end = std::next( face_verts_begin, num_face_verts );
        auto face_verts = std::vector<id_t>( face_verts_begin, std::next( face_verts_begin, num_face_verts ) );
        
        // get the edges of the face
        size_t num_face_edges = 0;
        auto face_edges_begin = from_domain_conn[2][1]->get_entities( face_id, num_face_edges );
        auto face_edges_end = std::next( face_edges_begin, num_face_edges );

        // get the cells of this face
        auto face_cells = from_domain_conn[2][3]->get_entities( face_id );
        // reverse the list of vertices if this is backwards
        if ( face_cells[0] != cell ) 
          std::reverse( face_verts.begin(), face_verts.end() );

        // a lambda function to locate an edge connected to two points
        auto _find_edge = [&]( const auto & pa, const auto & pb  ) 
        {
          // locate the edge with the two vertices
          auto edge = std::find_if( 
            face_edges_begin, face_edges_end, 
            [&]( const auto & e ) 
            { 
              size_t num_verts = 0;
              auto verts = from_domain_conn[1][0]->get_entities( e.entity(), num_verts );
              assert( num_verts == 2 && "should be two vertices per edge" );
              return ( (verts[0] == pa && verts[1] == pb) || 
                       (verts[0] == pb && verts[1] == pa) );
            } 
          );
          // make sure we found an edge
          assert( edge != face_edges_end );
          // return the edge
          return edge;
        };
          
        // loop over each edge (pair of vertices of the face)
        // there is a wedge for each vertex -> edge -> face combination
        auto p1 = std::prev( face_verts.end() );
        auto p0 = std::prev( p1 );
        auto edge0 = _find_edge( *p0, *p1 );
        for ( auto p2=face_verts.begin(); p2!=face_verts.end(); p2++ ) {
          // get the next edge
          auto edge1 = _find_edge( *p1, *p2 );
          // get the corner of this point, but also associated with this cell
          size_t num_point_corners = 0, num_cell_corners = 0;
          auto point_corners_begin = to_domain_conn[0][0]->get_entities( *p1, num_point_corners );
          auto cell_corners_begin  = to_domain_conn[3][0]->get_entities( cell_id, num_cell_corners );
          // copy the result for set intersection
          auto point_corners = std::vector<id_t>( point_corners_begin, std::next( point_corners_begin, num_point_corners ) );
          auto  cell_corners = std::vector<id_t>(  cell_corners_begin, std::next(  cell_corners_begin, num_cell_corners ) );
          // sort the lists for intersectinos
          std::sort( point_corners.begin(), point_corners.end() );
          std::sort(  cell_corners.begin(),  cell_corners.end() );
          // get the intersections of the sets
          std::vector<id_t> corners;
          corners.reserve(1);
          std::set_intersection( point_corners.begin(), point_corners.end(),
                                 cell_corners.begin(), cell_corners.end(),
                                 std::back_inserter(corners));
          assert( corners.size() == 1 );
          const auto & c1 = corners.front();
          // the first point, this is the right (even) one
          entities[i++] = *p1;
          entities[i++] = *edge0;
          entities[i++] = *face;
          entities[i++] = c1;
          // for the next point, this one is the left (odd) one
          entities[i++] = *p1;
          entities[i++] = *edge1;
          entities[i++] = *face;
          entities[i++] = c1;
          // update the iterators
          p0 = p1;
          p1 = p2;
          edge0 = edge1;
          // increment the wedge counter
          num_wedges += 2;
        }
      } // foreach vertex

      return std::vector<id_t>( num_wedges, 4 );
    }
      //------------------------------------------------------------------------
      // failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch
  } // create_bound_entities


  //============================================================================
  // Private Data
  //============================================================================
  
private:

  //! the list of directions for each face
  std::vector<char> face_direcions_;


};



} // namespace
} // namespace
