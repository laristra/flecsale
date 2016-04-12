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
#include "../../geom/geometric_shapes.h"
#include "../../math/math.h"
#include "../../mesh/burton/burton_entity_types.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_triangle_t burton_entity_types.h
//!
//! \brief The burton_triangle_t type provides a derived instance of
//!   burton_cell_t for 2D triangle cells.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_triangle_t : public burton_cell_t<N>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_cell_t<N>;

  //! the mesh topology type
  using mesh_topology_base_t = typename base_t::mesh_topology_base_t;

  //! the mesh traits
  using mesh_traits_t = typename base_t::mesh_traits_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  // the id type
  using id_t = typename base_t::id_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_triangle_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  // use some base accessors
  using base_t::edges;
  using base_t::vertices;

  //! the centroid
  point_t centroid() const override
  {
    auto vs = base_t::vertices();
    return math::average( 
      vs[0]->coordinates(), vs[1]->coordinates(), vs[2]->coordinates() );
  }

  //! the area of the cell
  real_t area() const override
  {
    auto vs = vertices();
    auto u = vs[1]->coordinates() - vs[0]->coordinates();
    auto v = vs[2]->coordinates() - vs[0]->coordinates();
    auto cross = cross_product( u, v );
    return std::abs( cross ) / 2;
  }

  //! the minimum length in the cell
  real_t min_length() const override
  {
    auto vs = vertices();
    // check the edges first
    auto es = edges();
    auto eit = es.begin();
    auto min_length = eit->length();
    std::for_each( ++eit, es.end(), [&](auto && e) 
                   { 
                     min_length = std::min( e->length(), min_length );
                   });
    // return result
    return min_length;
  }


  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::geometric_shapes_t::triangle; };

  //! \brief create_entities function for burton_triangle_cell_t.
  inline std::vector<id_t> create_entities(
      size_t dim, id_t * e, id_t * v, size_t vertex_count ) override
  {
    assert( vertex_count == 3 );
    e[0] = v[0];
    e[1] = v[1];

    e[2] = v[1];
    e[3] = v[2];

    e[4] = v[2];
    e[5] = v[0];

    return {2, 2, 2};
  } // create_entities

  //! \brief create_bound_entities function for burton_triangle_cell_t.
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, id_t ** ent_ids, 
    size_t * ent_counts, id_t * c )  override
  {
    assert( ent_counts[0] == 3 );
    switch (dim) {
      // Corners
      // The right edge is always first
      case 1:
        // corner 0
        c[0] = ent_ids[0][0]; // vertex 0
        c[1] = ent_ids[1][0]; // edge 0, abuts vertex 0
        c[2] = ent_ids[1][2]; // edge 3, abuts vertex 0

        // corner 1
        c[3] = ent_ids[0][1]; // vertex 1
        c[4] = ent_ids[1][1]; // edge 1, abuts vertex 1
        c[5] = ent_ids[1][0]; // edge 0, abuts vertex 1

        // corner 2
        c[6] = ent_ids[0][2]; // vertex 2
        c[7] = ent_ids[1][2]; // edge 2, abuts vertex 2
        c[8] = ent_ids[1][1]; // edge 1, abuts vertex 2

        return {3, 3, 3};

      // Wedges
      case 2:

        // wedge 0
        c[0] = ent_ids[0][0]; // vertex 0
        c[1] = ent_ids[1][2]; // edge 2

        // wedge 1
        c[2] = ent_ids[0][0]; // vertex 0
        c[3] = ent_ids[1][0]; // edge 0

        // wedge 2
        c[4] = ent_ids[0][1]; // vertex 1
        c[5] = ent_ids[1][0]; // edge 0

        // wedge 3
        c[6] = ent_ids[0][1]; // vertex 1
        c[7] = ent_ids[1][1]; // edge 1

        // wedge 4
        c[8] = ent_ids[0][2]; // vertex 2
        c[9] = ent_ids[1][1]; // edge 1

        // wedge 5
        c[10] = ent_ids[0][2]; // vertex 2
        c[11] = ent_ids[1][2]; // edge 2

        return {2, 2, 2, 2, 2, 2};

      default:
        assert(false && "Unknown bound entity type");
    } // switch
  } // create_bound_entities


}; // class burton_triangle_cell_t

} // namespace
} // namespace
