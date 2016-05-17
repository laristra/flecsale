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
#include "ale/mesh/burton/burton_mesh_topology.h"
#include "ale/mesh/burton/burton_polyhedron.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
// 3D polyhedron
////////////////////////////////////////////////////////////////////////////////

//! the centroid
burton_polyhedron_t::point_t burton_polyhedron_t::centroid() const
{
  auto coords = coordinates();
  return math::average( coords );
}

//! the centroid
burton_polyhedron_t::point_t burton_polyhedron_t::midpoint() const
{
  auto coords = coordinates();
  return math::average( coords );
}


//! the area of the cell
burton_polyhedron_t::real_t burton_polyhedron_t::volume() const
{
  auto coords = coordinates();
  return 0;
}


} // namespace
} // namespace
