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
#include "ale/geom/shapes/tetrahedron.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_tetrahedron_t burton_entity_types.h
//!
//! \brief The burton_tetrahedron_t type provides a derived instance of
//!   burton_cell_t for 2D tetrahedron cells.
////////////////////////////////////////////////////////////////////////////////
class burton_tetrahedron_t : public burton_element_t<3,3>
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
  burton_tetrahedron_t(mesh_topology_base_t & mesh) : base_t(mesh)
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
  { return geom::tetrahedron::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_tetrahedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_entities(
    size_t dim, const id_t & cell,
    connectivity_t*  (&conn)[num_dimensions+1][num_dimensions], 
    id_t * entities ) override
  {

    auto cell_id = cell.entity();
    size_t num_cell_verts = 0;
    auto v = conn[3][0]->get_entities( cell_id, num_cell_verts );

    assert( num_cell_verts == 4 );

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
      entities[5] = v[0];

      // top
      entities[6] = v[0];
      entities[7] = v[3];

      entities[8] = v[1];
      entities[9] = v[3];

      entities[10] = v[2];
      entities[11] = v[3];

      return {2, 2, 2, 2, 2, 2};

      //------------------------------------------------------------------------
      // Faces
    case (2): 

      entities[0]  = v[0];
      entities[1]  = v[1];
      entities[2]  = v[3];

      entities[3]  = v[1];
      entities[4]  = v[2];
      entities[5]  = v[3];

      entities[6]  = v[2];
      entities[7]  = v[0];
      entities[8]  = v[3];

      entities[9]  = v[0];
      entities[10] = v[2];
      entities[11] = v[1];

      return {3, 3, 3, 3};
      
      //------------------------------------------------------------------------
      // Failure      
    default: 
      raise_runtime_error("Unknown entity type");

    } // switch


  } // create_entities

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities function for burton_tetrahedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    connectivity_t*  (&from_domain_conn)[num_dimensions+1][num_dimensions+1], 
    connectivity_t*  (&  to_domain_conn)[num_dimensions+1][num_dimensions+1], 
    id_t * c )  override
  {
    size_t num_vertices = 0, num_edges = 0, num_faces = 0;
    auto verts = from_domain_conn[3][0]->get_entities( cell_id.entity(), num_vertices );
    auto edges = from_domain_conn[3][1]->get_entities( cell_id.entity(), num_edges );
    auto faces = from_domain_conn[3][2]->get_entities( cell_id.entity(), num_faces );

    assert( num_vertices == 4 );

    size_t i = 0;

    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
    case 0:

      // corner 0
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] =   edges[2]; // edge 2, abuts vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] =   faces[0]; // face 0, abuts vertex 0
      c[i++] =   faces[2]; // face 2, abuts vertex 0
      c[i++] =   faces[3]; // face 3, abuts vertex 0
      // corner 1
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] =   edges[4]; // edge 4, abuts vertex 1
      c[i++] =   faces[1]; // face 1, abuts vertex 1
      c[i++] =   faces[0]; // face 0, abuts vertex 1
      c[i++] =   faces[3]; // face 3, abuts vertex 1
      // corner 2
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] =   edges[5]; // edge 5, abuts vertex 2
      c[i++] =   faces[2]; // face 2, abuts vertex 2
      c[i++] =   faces[1]; // face 1, abuts vertex 2
      c[i++] =   faces[3]; // face 3, abuts vertex 2
      // corner 3
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] =   edges[5]; // edge 5, abuts vertex 3
      c[i++] =   edges[4]; // edge 4, abuts vertex 3
      c[i++] =   faces[0]; // face 0, abuts vertex 3
      c[i++] =   faces[1]; // face 1, abuts vertex 3
      c[i++] =   faces[2]; // face 2, abuts vertex 3

      return std::vector<id_t>(4, 7);

      //------------------------------------------------------------------------
      // Wedges
    case 1: {
      
      size_t num_corners = 0;
      auto corners = to_domain_conn[3][0]->get_entities( cell_id.entity(), num_corners );

      // corner 0
      // wedge 0
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] =   faces[0]; // face 0, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 1
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] =   faces[0]; // face 0, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 2
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[2]; // edge 2, abuts vertex 0
      c[i++] =   faces[2]; // face 2, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 3
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] =   faces[2]; // face 2, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 4
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] =   faces[3]; // face 3, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 5
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[2]; // edge 2, abuts vertex 0
      c[i++] =   faces[3]; // face 3, abuts vertex 0
      c[i++] = corners[0]; // corner 0

      // corner 1
      // wedge 0
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[4]; // edge 4, abuts vertex 1
      c[i++] =   faces[1]; // face 1, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 1
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] =   faces[1]; // face 1, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 2
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] =   faces[0]; // face 0, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 3
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[4]; // edge 4, abuts vertex 1
      c[i++] =   faces[0]; // face 0, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 4
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] =   faces[3]; // face 3, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 5
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] =   faces[3]; // face 3, abuts vertex 1
      c[i++] = corners[1]; // corner 1

      // corner 2
      // wedge 0
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[5]; // edge 5, abuts vertex 2
      c[i++] =   faces[2]; // face 2, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 1
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] =   faces[2]; // face 2, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 2
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] =   faces[1]; // face 1, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 3
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[5]; // edge 5, abuts vertex 2
      c[i++] =   faces[1]; // face 1, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 4
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] =   faces[3]; // face 3, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 5
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] =   faces[3]; // face 3, abuts vertex 2
      c[i++] = corners[2]; // corner 2


      // corner 3
      // wedge 0
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[4]; // edge 4, abuts vertex 3
      c[i++] =   faces[0]; // face 0, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 1
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] =   faces[0]; // face 0, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 2
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[5]; // edge 5, abuts vertex 3
      c[i++] =   faces[1]; // face 1, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 3
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[4]; // edge 4, abuts vertex 3
      c[i++] =   faces[1]; // face 1, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 4
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] =   faces[2]; // face 2, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 5
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[5]; // edge 5, abuts vertex 3
      c[i++] =   faces[2]; // face 2, abuts vertex 3
      c[i++] = corners[3]; // corner 3

      return std::vector<id_t>(24, 4);
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
