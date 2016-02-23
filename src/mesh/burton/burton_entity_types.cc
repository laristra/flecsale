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

//! user includes
#include "flecsi/mesh/mesh_topology.h"

#include "ale/math/operators.h"
#include "ale/mesh/burton/burton_types.h"

namespace ale {
namespace mesh {

//!  the real type
using real_t = burton_mesh_traits_t::real_t;

//!  the vector type
using vector_t = burton_mesh_traits_t::vector_t;

//!  the point type
using point_t = burton_mesh_traits_t::point_t;

/*----------------------------------------------------------------------------*
 * burton_edge_t
 *----------------------------------------------------------------------------*/

point_t burton_edge_t::midpoint() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    return point_t{0.5*(vs[0]->coordinates() + vs[1]->coordinates())};
  } // burton_edge_t::midpoint

real_t burton_edge_t::length() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    auto & a = vs[0]->coordinates();
    auto & b = vs[1]->coordinates();
    
    using math::sqr;

    return std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
  } // burton_edge_t::length

vector_t burton_edge_t::normal() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    auto & a = vs[0]->coordinates();
    auto & b = vs[1]->coordinates();

    return { a[1] - b[1], b[0] - a[0] };
  } // burton_edge_t::normal

/*----------------------------------------------------------------------------*
 * burton_quadrilateral_cell_t
 *----------------------------------------------------------------------------*/

point_t burton_quadrilateral_cell_t::centroid() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    auto tmp = vs[0]->coordinates();
    tmp += vs[1]->coordinates();
    tmp += vs[2]->coordinates();
    tmp += vs[3]->coordinates();
    tmp /= 4;

    return tmp;
  } // burton_quadrilateral_cell_t::centroid

real_t burton_quadrilateral_cell_t::area() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    auto p = vs[0]->coordinates() - vs[2]->coordinates();
    auto q = vs[1]->coordinates() - vs[3]->coordinates();
    
    using math::cross_product;

    auto det = cross_product( p, q );

    return 0.5 * std::abs( det );
  } // burton_quadrilateral_cell_t::area

} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
