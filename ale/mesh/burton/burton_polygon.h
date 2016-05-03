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

  //! the area of the cell
  real_t area() const override;

  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::polygon<num_dimensions>::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_entities(
      size_t dim, id_t * e, id_t * v, size_t vertex_count) override
  {
    assert( dim == 1 );
    auto vp = v[vertex_count - 1];
    for ( auto i=0, ind=0; i<vertex_count; i++ ) {
      auto vn = v[i];
      e[ ind++ ] = vp;
      e[ ind++ ] = vn;
      vp = vn;
    }
    return std::vector<id_t>(vertex_count, 2);
  }

  //----------------------------------------------------------------------------
  //!  \brief create_bound_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, id_t ** ent_ids, 
    size_t * ent_counts, id_t * c ) override
  {
    auto vertex_count = ent_counts[0];

    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      // The right edge is always first
    case 1: {

      auto vp = vertex_count - 1;
      for ( auto i=0, ind=0; i<vertex_count; i++ ) {
        auto vn = i;
        c[ ind++ ] = ent_ids[0][vn]; // vertex 0
        c[ ind++ ] = ent_ids[1][vn]; // edge 0, abuts vertex 0
        c[ ind++ ] = ent_ids[1][vp]; // edge 3, abuts vertex 0
        vp = vn;
      }
      return std::vector<id_t>(vertex_count, 3);
    }

      //------------------------------------------------------------------------
      // Wedges
    case 2: {
      
      auto vp = vertex_count - 1;
      for ( auto i=0, ind=0; i<vertex_count; i++ ) {
        auto vn = i;
        // wedge 0
        c[ ind++ ] = ent_ids[0][vn]; // vertex 0
        c[ ind++ ] = ent_ids[1][vp]; // edge 3
        // wedge 1
        c[ ind++ ] = ent_ids[0][vn]; // vertex 0
        c[ ind++ ] = ent_ids[1][vn]; // edge 0
        vp = vn;
      }
      return std::vector<id_t>(2*vertex_count, 2);
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
  inline std::vector<id_t> create_entities(
      size_t dim, id_t * e, id_t * v, size_t vertex_count) override
  {
    assert( dim == 1 );
    auto vp = v[vertex_count - 1];
    for ( auto i=0, ind=0; i<vertex_count; i++ ) {
      auto vn = v[i];
      e[ ind++ ] = vp;
      e[ ind++ ] = vn;
      vp = vn;
    }
    return std::vector<id_t>(vertex_count, 2);
  }

  //----------------------------------------------------------------------------
  //!  \brief create_bound_entities function for burton_polygonal_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, id_t ** ent_ids, 
    size_t * ent_counts, id_t * c ) override
  {
    auto vertex_count = ent_counts[0];

    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      // The right edge is always first
    case 1: {

      auto vp = vertex_count - 1;
      for ( auto i=0, ind=0; i<vertex_count; i++ ) {
        auto vn = i;
        c[ ind++ ] = ent_ids[0][vn]; // vertex 0
        c[ ind++ ] = ent_ids[1][vn]; // edge 0, abuts vertex 0
        c[ ind++ ] = ent_ids[1][vp]; // edge 3, abuts vertex 0
        vp = vn;
      }
      return std::vector<id_t>(vertex_count, 3);
    }

      //------------------------------------------------------------------------
      // Wedges
    case 2: {
      
      auto vp = vertex_count - 1;
      for ( auto i=0, ind=0; i<vertex_count; i++ ) {
        auto vn = i;
        // wedge 0
        c[ ind++ ] = ent_ids[0][vn]; // vertex 0
        c[ ind++ ] = ent_ids[1][vp]; // edge 3
        // wedge 1
        c[ ind++ ] = ent_ids[0][vn]; // vertex 0
        c[ ind++ ] = ent_ids[1][vn]; // edge 0
        vp = vn;
      }
      return std::vector<id_t>(2*vertex_count, 2);
    }
      //------------------------------------------------------------------------
      // Failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch

  } // create_bound_entities


};


} // namespace
} // namespace
