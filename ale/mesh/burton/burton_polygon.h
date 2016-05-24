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
#include "ale/geom/shapes/polygon.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_polygon_t burton_entity_types.h
//!
//! \brief The burton_polygon_t type provides a derived instance of
//!   burton_cell_t for 2D polygon cells.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_polygon_t {};

////////////////////////////////////////////////////////////////////////////////
//! \class burton_polygon_t burton_entity_types.h
//!
//! \brief The burton_polygon_t type provides a derived instance of
//!   burton_cell_t for 2D polygon cells.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_polygon_t<2> : public burton_element_t<2,2>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<2,2>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_polygon_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  point_t centroid() const override;

  //! the midpoint
  point_t midpoint() const override;

  //! the area of the cell
  real_t area() const override;

  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::polygon<num_dimensions>::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_entities(
    size_t dim, const id_t & cell,
    connectivity_t*  (&conn)[num_dimensions+1][num_dimensions], 
    id_t * entities ) override
  {
    assert( dim == 1 );

    auto cell_id = cell.entity();
    size_t num_cell_verts = 0;
    auto v = conn[2][0]->get_entities( cell_id, num_cell_verts );

    size_t ind=0;
    for ( auto i=0; i<num_cell_verts-1; i++ ) {
      auto vp = v[i];
      auto vn = v[i+1];
      entities[ ind++ ] = vp;
      entities[ ind++ ] = vn;
    }
    entities[ ind++ ] = v[ num_cell_verts-1 ];
    entities[ ind++ ] = v[ 0 ];

    return std::vector<id_t>(num_cell_verts, 2);
  }

  //----------------------------------------------------------------------------
  //!  \brief create_bound_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    connectivity_t*  (&from_domain_conn)[num_dimensions+1][num_dimensions+1], 
    connectivity_t*  (&  to_domain_conn)[num_dimensions+1][num_dimensions+1], 
    id_t * c )  override
  {
    size_t num_vertices = 0, num_edges = 0, num_corners = 0;
    auto verts = from_domain_conn[2][0]->get_entities( cell_id.entity(), num_vertices );
    auto edges = from_domain_conn[2][1]->get_entities( cell_id.entity(), num_edges );

    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      // The right edge is always first
    case 0: {

      auto vp = num_vertices - 1;
      for ( auto vn=0, ind=0; vn<num_vertices; vn++ ) {
        c[ ind++ ] = verts[vn]; // vertex 0
        c[ ind++ ] = edges[vn]; // edge 0, abuts vertex 0
        c[ ind++ ] = edges[vp]; // edge 3, abuts vertex 0
        vp = vn;
      }
      return std::vector<id_t>(num_vertices, 3);
    }
      //------------------------------------------------------------------------
      // wedges
      // The right wedge is always first
    case 1: {

      size_t num_corners = 0;
      auto corners = to_domain_conn[2][0]->get_entities( cell_id.entity(), num_corners );

      auto vp = num_vertices - 1;
      for ( auto vn=0, ind=0; vn<num_vertices; vn++ ) {
        // wedge 0
        c[ ind++ ] =   verts[vn]; // vertex 0
        c[ ind++ ] =   edges[vn]; // edge 0, abuts vertex 0
        c[ ind++ ] = corners[vn]; // corner 0
        // wedge 1
        c[ ind++ ] =   verts[vn]; // vertex 0
        c[ ind++ ] =   edges[vp]; // edge 3, abuts vertex 0
        c[ ind++ ] = corners[vn]; // corner 0
        vp = vn;
      }
      return std::vector<id_t>(2*num_vertices, 3);
    }
      //------------------------------------------------------------------------
      // Failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch

  } // create_bound_entities


};


////////////////////////////////////////////////////////////////////////////////
//! \class burton_polygon_t burton_entity_types.h
//!
//! \brief The burton_polygon_t type provides a derived instance of
//!   burton_cell_t for 2D polygon cells.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_polygon_t<3> : public burton_element_t<3,2>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<3,2>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_polygon_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  point_t centroid() const override;

  //! the midpoint
  point_t midpoint() const override;

  //! the normal
  vector_t normal() const override;

  //! the area of the cell
  real_t area() const override;

  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::polygon<num_dimensions>::shape; };

};


} // namespace
} // namespace
