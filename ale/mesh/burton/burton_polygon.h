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
#include "ale/geom/area.h"
#include "ale/geom/centroid.h"
#include "ale/geom/geometric_shapes.h"
#include "ale/mesh/burton/burton_entity_types.h"

namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_polygonal_t burton_entity_types.h
//!
//! \brief The burton_polygonal_t type provides a derived instance of
//!   burton_cell_t for 2D polygonal cells.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_polygon_t : public burton_element_t<N,2>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<N,2>;

  //! the flecsi mesh topology type
  using typename base_t::mesh_topology_base_t;

  //! the mesh traits
  using typename base_t::mesh_traits_t;

  //! Type containing coordinates of the vertex.
  using typename base_t::point_t;

  //! Type of floating point.
  using typename base_t::real_t;

  // the id type
  using typename base_t::id_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_polygon_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  // use some base accessors
  using base_t::coordinates;

  //! the centroid
  point_t centroid() const override
  {
    auto coords = coordinates();
    return geom::centroid( coords );
  }

  //! the area of the cell
  real_t area() const override
  {
    auto coords = coordinates();
    return geom::area( coords );
  }

  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::geometric_shapes_t::polygon; };

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
