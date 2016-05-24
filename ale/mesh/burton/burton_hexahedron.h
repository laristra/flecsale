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
#include "ale/geom/shapes/hexahedron.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_hexahedron_t burton_entity_types.h
//!
//! \brief The burton_hexahedron_t type provides a derived instance of
//!   burton_cell_t for 2D hexahedron cells.
////////////////////////////////////////////////////////////////////////////////
class burton_hexahedron_t : public burton_element_t<3,3>
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
  burton_hexahedron_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

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
  { return geom::hexahedron::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_hexahedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_entities(
    size_t dim, const id_t & cell,
    connectivity_t*  (&conn)[num_dimensions+1][num_dimensions], 
    id_t * entities ) override
  {

    auto cell_id = cell.entity();
    size_t num_cell_verts = 0;
    auto v = conn[3][0]->get_entities( cell_id, num_cell_verts );

    assert( num_cell_verts == 8 );
    
    switch (dim) {

      //------------------------------------------------------------------------
      // Edges
    case (1): 
      // bottom
      entities[0] = v[0];
      entities[1] = v[1];

      entities[2] = v[1];
      entities[3] = v[2];

      entities[4] = v[2];
      entities[5] = v[3];

      entities[6] = v[3];
      entities[7] = v[0];

      // top
      entities[8] = v[4];
      entities[9] = v[5];

      entities[10] = v[5];
      entities[11] = v[6];

      entities[12] = v[6];
      entities[13] = v[7];

      entities[14] = v[7];
      entities[15] = v[4];

      // vertical ones
      entities[16] = v[0];
      entities[17] = v[4];

      entities[18] = v[1];
      entities[19] = v[5];

      entities[20] = v[2];
      entities[21] = v[6];

      entities[22] = v[3];
      entities[23] = v[7];

      return {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

      //------------------------------------------------------------------------
      // Faces
    case (2): 
      // bottom
      entities[0]  = v[3];
      entities[1]  = v[2];
      entities[2]  = v[1];
      entities[3]  = v[0];
      // top
      entities[4]  = v[5];
      entities[5]  = v[6];
      entities[6]  = v[7];
      entities[7]  = v[4];
      // front
      entities[8]  = v[1];
      entities[9]  = v[5];
      entities[10] = v[4];
      entities[11] = v[0];
      // right
      entities[12] = v[2];
      entities[13] = v[6];
      entities[14] = v[5];
      entities[15] = v[1];
      // back
      entities[16] = v[3];
      entities[17] = v[7];
      entities[18] = v[6];
      entities[19] = v[2];
      // left
      entities[20] = v[0];
      entities[21] = v[4];
      entities[22] = v[7];
      entities[23] = v[3];

      return {4, 4, 4, 4, 4, 4};
      
      //------------------------------------------------------------------------
      // Failure      
    default: 
      raise_runtime_error("Unknown entity type");

    } // switch


  } // create_entities

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities function for burton_hexahedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    connectivity_t*  (&from_domain_conn)[num_dimensions+1][num_dimensions+1], 
    connectivity_t*  (&  to_domain_conn)[num_dimensions+1][num_dimensions+1], 
    id_t * c ) override
  {
    size_t num_vertices = 0, num_edges = 0, num_faces = 0;
    auto verts = from_domain_conn[3][0]->get_entities( cell_id.entity(), num_vertices );
    auto edges = from_domain_conn[3][1]->get_entities( cell_id.entity(), num_edges );
    auto faces = from_domain_conn[3][2]->get_entities( cell_id.entity(), num_faces );
    
    assert( num_vertices == 8 );

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
    case 0:
      // corner 0
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] =   edges[8]; // edge 8, abuts vertex 0
      c[i++] =   faces[2]; // face 2, abuts vertex 0
      c[i++] =   faces[5]; // face 5, abuts vertex 0
      c[i++] =   faces[0]; // face 0, abuts vertex 0

      // corner 1
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] =   edges[9]; // edge 9, abuts vertex 1
      c[i++] =   faces[3]; // face 3, abuts vertex 1
      c[i++] =   faces[2]; // face 2, abuts vertex 1
      c[i++] =   faces[0]; // face 0, abuts vertex 1

      // corner 2
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] =   edges[10]; // edge 10, abuts vertex 2
      c[i++] =   faces[4]; // face 4, abuts vertex 2
      c[i++] =   faces[3]; // face 3, abuts vertex 2
      c[i++] =   faces[0]; // face 0, abuts vertex 2

      // corner 3
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] =   edges[2]; // edge 2, abuts vertex 3
      c[i++] =   edges[11]; // edge 11, abuts vertex 3
      c[i++] =   faces[5]; // face 5, abuts vertex 3
      c[i++] =   faces[4]; // face 4, abuts vertex 3
      c[i++] =   faces[0]; // face 0, abuts vertex 3

      // corner 4
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[7]; // edge 7, abuts vertex 4
      c[i++] =   edges[4]; // edge 4, abuts vertex 4
      c[i++] =   edges[8]; // edge 8, abuts vertex 4
      c[i++] =   faces[5]; // face 5, abuts vertex 4
      c[i++] =   faces[2]; // face 2, abuts vertex 4
      c[i++] =   faces[1]; // face 1, abuts vertex 4

      // corner 5
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[4]; // edge 4, abuts vertex 5
      c[i++] =   edges[5]; // edge 5, abuts vertex 5
      c[i++] =   edges[9]; // edge 9, abuts vertex 5
      c[i++] =   faces[2]; // face 2, abuts vertex 5
      c[i++] =   faces[3]; // face 3, abuts vertex 5
      c[i++] =   faces[1]; // face 1, abuts vertex 5

      // corner 6
      c[i++] =   verts[6]; // vertex 6
      c[i++] =   edges[5]; // edge 5, abuts vertex 6
      c[i++] =   edges[6]; // edge 6, abuts vertex 6
      c[i++] =   edges[10]; // edge 10, abuts vertex 6
      c[i++] =   faces[3]; // face 3, abuts vertex 6
      c[i++] =   faces[4]; // face 4, abuts vertex 6
      c[i++] =   faces[1]; // face 1, abuts vertex 6

      // corner 7
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[6]; // edge 6, abuts vertex 7
      c[i++] =   edges[7]; // edge 7, abuts vertex 7
      c[i++] =   edges[11]; // edge 11, abuts vertex 7
      c[i++] =   faces[4]; // face 4, abuts vertex 7
      c[i++] =   faces[5]; // face 5, abuts vertex 7
      c[i++] =   faces[1]; // face 1, abuts vertex 7

      return std::vector<id_t>(8, 7);

      //------------------------------------------------------------------------
      // Wedges
      //
      // Each corner has 6 vertices.  There is an even/odd ordering so we 
      // no which way to compute normals.  So edges are defined in pairs
      //
    case 1: {
  
      size_t num_corners = 0;
      auto corners = to_domain_conn[3][0]->get_entities( cell_id.entity(), num_corners );

      // corner 0
      // wedge 0
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[8]; // edge 8, abuts vertex 0
      c[i++] =   faces[2]; // face 2, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 1
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] =   faces[2]; // face 2, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 2
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] =   faces[5]; // face 5, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 3
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[8]; // edge 8, abuts vertex 0
      c[i++] =   faces[5]; // face 5, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 4
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] =   faces[0]; // face 0, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 5
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] =   faces[0]; // face 0, abuts vertex 0
      c[i++] = corners[0]; // corner 0

      // corner 1
      // wedge 0
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[9]; // edge 9, abuts vertex 1
      c[i++] =   faces[3]; // face 3, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 1
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] =   faces[3]; // face 3, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 2
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] =   faces[2]; // face 2, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 3
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[9]; // edge 9, abuts vertex 1
      c[i++] =   faces[2]; // face 2, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 4
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] =   faces[0]; // face 0, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 5
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] =   faces[0]; // face 0, abuts vertex 1
      c[i++] = corners[1]; // corner 1

      // corner 2
      // wedge 0
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[10]; // edge 10, abuts vertex 2
      c[i++] =   faces[4]; // face 4, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 1
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] =   faces[4]; // face 4, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 2
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] =   faces[3]; // face 3, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 3
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[10]; // edge 10, abuts vertex 2
      c[i++] =   faces[3]; // face 3, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 4
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] =   faces[0]; // face 0, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 5
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] =   faces[0]; // face 0, abuts vertex 2
      c[i++] = corners[2]; // corner 2

      // corner 3
      // wedge 0
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[11]; // edge 11, abuts vertex 3
      c[i++] =   faces[5]; // face 5, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 1
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] =   faces[5]; // face 5, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 2
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[2]; // edge 2, abuts vertex 3
      c[i++] =   faces[4]; // face 4, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 3
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[11]; // edge 11, abuts vertex 3
      c[i++] =   faces[4]; // face 4, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 4
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] =   faces[0]; // face 0, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 5
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[2]; // edge 2, abuts vertex 3
      c[i++] =   faces[0]; // face 0, abuts vertex 3
      c[i++] = corners[3]; // corner 3

      // corner 4
      // wedge 0
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[8]; // edge 8, abuts vertex 4
      c[i++] =   faces[5]; // face 5, abuts vertex 4
      c[i++] = corners[4]; // corner 4
      // wedge 1
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[7]; // edge 7, abuts vertex 4
      c[i++] =   faces[5]; // face 5, abuts vertex 4
      c[i++] = corners[4]; // corner 4
      // wedge 2
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[4]; // edge 4, abuts vertex 4
      c[i++] =   faces[2]; // face 2, abuts vertex 4
      c[i++] = corners[4]; // corner 4
      // wedge 3
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[8]; // edge 8, abuts vertex 4
      c[i++] =   faces[2]; // face 2, abuts vertex 4
      c[i++] = corners[4]; // corner 4
      // wedge 4
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[7]; // edge 7, abuts vertex 4
      c[i++] =   faces[1]; // face 1, abuts vertex 4
      c[i++] = corners[4]; // corner 4
      // wedge 5
      c[i++] =   verts[4]; // vertex 4
      c[i++] =   edges[4]; // edge 4, abuts vertex 4
      c[i++] =   faces[1]; // face 1, abuts vertex 4
      c[i++] = corners[4]; // corner 4

      // corner 5
      // wedge 0
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[9]; // edge 9, abuts vertex 5
      c[i++] =   faces[2]; // face 2, abuts vertex 5
      c[i++] = corners[5]; // corner 5
      // wedge 1
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[4]; // edge 4, abuts vertex 5
      c[i++] =   faces[2]; // face 2, abuts vertex 5
      c[i++] = corners[5]; // corner 5
      // wedge 2
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[5]; // edge 5, abuts vertex 5
      c[i++] =   faces[3]; // face 3, abuts vertex 5
      c[i++] = corners[5]; // corner 5
      // wedge 3
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[9]; // edge 9, abuts vertex 5
      c[i++] =   faces[3]; // face 3, abuts vertex 5
      c[i++] = corners[5]; // corner 5
      // wedge 4
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[4]; // edge 4, abuts vertex 5
      c[i++] =   faces[1]; // face 1, abuts vertex 5
      c[i++] = corners[5]; // corner 5
      // wedge 5
      c[i++] =   verts[5]; // vertex 5
      c[i++] =   edges[5]; // edge 5, abuts vertex 5
      c[i++] =   faces[1]; // face 1, abuts vertex 5
      c[i++] = corners[5]; // corner 5

      // corner 6
      // wedge 0
      c[i++] =   verts[6]; // vertex 6
      c[i++] =   edges[10]; // edge 10, abuts vertex 6
      c[i++] =   faces[3]; // face 3, abuts vertex 6
      c[i++] = corners[6]; // corner 6
      // wedge 1
      c[i++] =   verts[6]; // vertex 6
      c[i++] =   edges[5]; // edge 5, abuts vertex 6
      c[i++] =   faces[3]; // face 3, abuts vertex 6
      c[i++] = corners[6]; // corner 6
      // wedge 2
      c[i++] =   verts[6]; // vertex 6 
      c[i++] =   edges[6]; // edge 6, abuts vertex 6
      c[i++] =   faces[4]; // face 4, abuts vertex 6
      c[i++] = corners[6]; // corner 6
      // wedge 3
      c[i++] =   verts[6]; // vertex 6
      c[i++] =   edges[10]; // edge 10, abuts vertex 6
      c[i++] =   faces[4]; // face 4, abuts vertex 6
      c[i++] = corners[6]; // corner 6
      // wedge 4
      c[i++] =   verts[6]; // vertex 6
      c[i++] =   edges[5]; // edge 5, abuts vertex 6
      c[i++] =   faces[1]; // face 1, abuts vertex 6
      c[i++] = corners[6]; // corner 6
      // wedge 5
      c[i++] =   verts[6]; // vertex 6
      c[i++] =   edges[6]; // edge 6, abuts vertex 6
      c[i++] =   faces[1]; // face 1, abuts vertex 6
      c[i++] = corners[6]; // corner 6

      // corner 7
      // wedge 0
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[11]; // edge 11, abuts vertex 7
      c[i++] =   faces[4]; // face 4, abuts vertex 7
      c[i++] = corners[7]; // corner 7
      // wedge 1
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[6]; // edge 6, abuts vertex 7
      c[i++] =   faces[4]; // face 4, abuts vertex 7
      c[i++] = corners[7]; // corner 7
      // wedge 2
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[7]; // edge 7, abuts vertex 7
      c[i++] =   faces[5]; // face 5, abuts vertex 7
      c[i++] = corners[7]; // corner 7
      // wedge 3
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[11]; // edge 11, abuts vertex 7
      c[i++] =   faces[5]; // face 5, abuts vertex 7
      c[i++] = corners[7]; // corner 7
      // wedge 4
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[6]; // edge 6, abuts vertex 7
      c[i++] =   faces[1]; // face 1, abuts vertex 7
      c[i++] = corners[7]; // corner 7
      // wedge 5
      c[i++] =   verts[7]; // vertex 7
      c[i++] =   edges[7]; // edge 7, abuts vertex 7
      c[i++] =   faces[1]; // face 1, abuts vertex 7
      c[i++] = corners[7]; // corner 7

      return std::vector<id_t>(48, 4);
    }
      //------------------------------------------------------------------------
      // failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch
  } // create_bound_entities

};



} // namespace
} // namespace
