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

#include "ale/geom/centroid.h"
#include "ale/math/operators.h"
#include "ale/mesh/burton/burton_types.h"

namespace ale {
namespace mesh {

//!  the point type
using point_t = burton_cell_t::point_t;


/*----------------------------------------------------------------------------*
 * burton_edge_t
 *----------------------------------------------------------------------------*/

point_t burton_edge_t::midpoint() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    return math::average( vs[0]->coordinates(), vs[1]->coordinates() );
  } // burton_edge_t::midpoint

/*----------------------------------------------------------------------------*
 * burton_quadrilateral_cell_t
 *----------------------------------------------------------------------------*/

point_t burton_quadrilateral_cell_t::centroid() const
  {
    auto & mesh = static_cast<mesh_topology_t<burton_mesh_types_t> &>(mesh_);
    auto vs = mesh.entities<0,0>(this).to_vec();

    return geom::centroid( vs[0]->coordinates(), vs[1]->coordinates(), 
                           vs[2]->coordinates(), vs[3]->coordinates() );
  } // burton_quadrilateral_cell_t::centroid

} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
