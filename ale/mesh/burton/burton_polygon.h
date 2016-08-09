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
  inline std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities ) override
  {
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, vertex_t::dimension );
    auto num_cell_verts = v.size();

    size_t ind=0;
    for ( auto i=0; i<num_cell_verts-1; i++ ) {
      auto vp = v[i];
      auto vn = v[i+1];
      entities[ ind++ ] = vp;
      entities[ ind++ ] = vn;
    }
    entities[ ind++ ] = v[ num_cell_verts-1 ];
    entities[ ind++ ] = v[ 0 ];

    return std::vector<size_t>(num_cell_verts, 2);
  }

  //----------------------------------------------------------------------------
  //!  \brief create_bound_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * c )  override
  {

    auto verts = primal_conn.get_entity_vec( cell, vertex_t::dimension );
    auto edges = primal_conn.get_entity_vec( cell,   edge_t::dimension );
    auto num_vertices = verts.size();

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
      return std::vector<size_t>(num_vertices, 3);
    }
      //------------------------------------------------------------------------
      // wedges
      // The right wedge is always first
    case 1: {

      auto corners = domain_conn.get_entity_vec( cell, corner_t::dimension );

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
      return std::vector<size_t>(2*num_vertices, 3);
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


  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities ) override
  {
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, vertex_t::dimension );
    auto num_cell_verts = v.size();

    size_t ind=0;
    for ( auto i=0; i<num_cell_verts-1; i++ ) {
      auto vp = v[i];
      auto vn = v[i+1];
      entities[ ind++ ] = vp;
      entities[ ind++ ] = vn;
    }
    entities[ ind++ ] = v[ num_cell_verts-1 ];
    entities[ ind++ ] = v[ 0 ];

    return std::vector<size_t>(num_cell_verts, 2);
  }

};


} // namespace
} // namespace
