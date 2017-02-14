/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines polyhedral elements for the burton_mesh_t class.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "geom/shapes/polyhedron.h"
#include "mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_polyhedron_t type provides a derived instance of
//!   burton_cell_t for 3D polyhedron cells.
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
  geom::shapes::geometric_shapes_t type() const override 
  { return geom::shapes::polyhedron<point_t>::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_polyhedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities )  override
  {        

    // you should only be coming in here to build edges
    assert( dim == 1 ); 

    // get the cell entities
    auto cell_faces = conn.get_entity_vec( cell, face_t::dimension );
    // make sure the faces exist
    assert( cell_faces.size() > 0 && "no cell faces yet" );
    
    // the list of edge pairs
    std::vector< std::pair< id_t, id_t > > cell_edges;

    // get the edges of each face
    for ( const auto & face : cell_faces ) { 
      // get all vertices in this face
      auto face_verts = conn.get_entity_vec( face, vertex_t::dimension );
      assert( face_verts.size() > 2 && "not enough vertices for a valid face" );
      // reserve space 
      cell_edges.reserve( cell_edges.size() + face_verts.size() );
      // add each edge pair to the list
      auto vp = std::prev( face_verts.end() );
      assert( vp != face_verts.begin() && "no vertices in this face" ); 
      for ( auto vn = face_verts.begin(); vn != face_verts.end(); vn++ ) {
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
        id_t a1 = a.first;
        id_t b1 = b.first;
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

    return std::vector<size_t>( num_edges, 2 );

      
  }

  //----------------------------------------------------------------------------
  // \brief create_bound_entities function for burton_polyhedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities )  override
  {

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

      std::vector<size_t> entity_count;

      // get the cell entities
      auto cell_verts = primal_conn.get_entity_vec( cell, vertex_t::dimension );
      auto cell_edges = primal_conn.get_entity_vec( cell, edge_t::dimension ).vec();
      auto cell_faces = primal_conn.get_entity_vec( cell, face_t::dimension ).vec();

      // sort the edges and faces for intersections later      
      std::sort( cell_edges.begin(), cell_edges.end() );
      std::sort( cell_faces.begin(), cell_faces.end() );

      // temparary lists
      std::vector<id_t> edges, faces;

      for ( const auto & vert : cell_verts ) { 
        // clear temporary lits
        edges.clear();
        faces.clear();
        // get all entities attached to this vertex
        auto vert_edges = primal_conn.get_entity_vec( vert, edge_t::dimension ).vec(); 
        auto vert_faces = primal_conn.get_entity_vec( vert, face_t::dimension ).vec(); 
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
        entities[i++] = vert;
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
      auto cell_faces   = primal_conn.get_entity_vec( cell,   face_t::dimension );
      auto cell_corners = domain_conn.get_entity_vec( cell, corner_t::dimension );
      
      // a counter for the number of wedges
      size_t num_wedges = 0;

      // loop over faces
      for ( const auto & face : cell_faces ) { 

        // get the vertices of the face
        auto face_verts = primal_conn.get_entity_vec( face, vertex_t::dimension );
        auto ccw_face_verts = std::vector<id_t>( face_verts.begin(), face_verts.end() );
        
        // get the edges of the face
        auto face_edges = primal_conn.get_entity_vec( face, edge_t::dimension );

        // get the cells of this face
        auto face_cells = primal_conn.get_entity_vec( face, cell_t::dimension );
        // reverse the list of vertices if this is backwards
        if ( face_cells[0] != cell ) 
          std::reverse( ccw_face_verts.begin(), ccw_face_verts.end() );

        // a lambda function to locate an edge connected to two points
        auto _find_edge = [&]( const auto & pa, const auto & pb  ) 
        {
          // locate the edge with the two vertices
          auto edge = std::find_if( 
            face_edges.begin(), face_edges.end(), 
            [&]( const auto & e ) 
            { 
              auto verts = primal_conn.get_entity_vec( e, vertex_t::dimension );
              assert( verts.size() == 2 && "should be two vertices per edge" );
              return ( (verts[0] == pa && verts[1] == pb) || 
                       (verts[0] == pb && verts[1] == pa) );
            } 
          );
          // make sure we found an edge
          assert( edge != face_edges.end() );
          // return the edge
          return edge;
        };
          
        // loop over each edge (pair of vertices of the face)
        // there is a wedge for each vertex -> edge -> face combination
        auto p1 = std::prev( ccw_face_verts.end() );
        auto p0 = std::prev( p1 );
        auto edge0 = _find_edge( *p0, *p1 );
        for ( auto p2=ccw_face_verts.begin(); p2!=ccw_face_verts.end(); p2++ ) {
          // get the next edge
          auto edge1 = _find_edge( *p1, *p2 );
          // get the corner of this point, but also associated with this cell
          auto point_corners = domain_conn.get_entity_vec(  *p1, corner_t::dimension ).vec();
          auto cell_corners  = domain_conn.get_entity_vec( cell, corner_t::dimension ).vec();
          // sort the lists for intersectinos
          std::sort( point_corners.begin(), point_corners.end() );
          std::sort( cell_corners.begin(),  cell_corners.end() );
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
          entities[i++] = face;
          entities[i++] = c1;
          // for the next point, this one is the left (odd) one
          entities[i++] = *p1;
          entities[i++] = *edge1;
          entities[i++] = face;
          entities[i++] = c1;
          // update the iterators
          p0 = p1;
          p1 = p2;
          edge0 = edge1;
          // increment the wedge counter
          num_wedges += 2;
        }
      } // foreach vertex

      return std::vector<size_t>( num_wedges, 4 );
    }
      //------------------------------------------------------------------------
      // failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch
  }


  //============================================================================
  // Private Data
  //============================================================================
  
private:

  //! the list of directions for each face
  std::vector<char> face_direcions_;


};



} // namespace
} // namespace
