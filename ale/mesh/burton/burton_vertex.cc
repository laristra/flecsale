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
#include "ale/mesh/burton/burton_vertex.h"
#include "ale/mesh/burton/burton_mesh_topology.h"


namespace ale {
namespace mesh {

// some type aliases
using burton_2d_vertex_t = burton_vertex_t<2>;
using burton_3d_vertex_t = burton_vertex_t<3>;


////////////////////////////////////////////////////////////////////////////////
// Get the coordinates at a vertex from the state handle.
////////////////////////////////////////////////////////////////////////////////
const burton_2d_vertex_t::point_t & burton_2d_vertex_t::coordinates() const
{
  using flecsi::mesh_entity_base_t;
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  const auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
    "coordinates", mesh->runtime_id() );
  return c[mesh_entity_base_t<num_domains>::template id<0>()];
} // coordinates

const burton_3d_vertex_t::point_t & burton_3d_vertex_t::coordinates() const
{
  using flecsi::mesh_entity_base_t;
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  const auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
    "coordinates", mesh->runtime_id() );
  return c[mesh_entity_base_t<num_domains>::template id<0>()];
} // coordinates

////////////////////////////////////////////////////////////////////////////////
// Get the coordinates at a vertex from the state handle.
////////////////////////////////////////////////////////////////////////////////
burton_2d_vertex_t::point_t & burton_2d_vertex_t::coordinates()
{
  using flecsi::mesh_entity_base_t;
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
    "coordinates", mesh->runtime_id() );
  return c[mesh_entity_base_t<num_domains>::template id<0>()];
} // coordinates

burton_3d_vertex_t::point_t & burton_3d_vertex_t::coordinates()
{
  using flecsi::mesh_entity_base_t;
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
    "coordinates", mesh->runtime_id() );
  return c[mesh_entity_base_t<num_domains>::template id<0>()];
} // coordinates

////////////////////////////////////////////////////////////////////////////////
// is this a boundary
////////////////////////////////////////////////////////////////////////////////
bool burton_2d_vertex_t::is_boundary() const
{
  using flecsi::mesh_entity_base_t;
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto flag =
    data_t::instance().template dense_accessor<bitfield_t, flecsi_internal>(
      "point_flags", mesh->runtime_id() );
  return flag[mesh_entity_base_t<num_domains>::template id<0>()].anybitset();
}

bool burton_3d_vertex_t::is_boundary() const
{
  using flecsi::mesh_entity_base_t;
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto flag =
    data_t::instance().template dense_accessor<bitfield_t, flecsi_internal>(
      "point_flags", mesh->runtime_id() );
  return flag[mesh_entity_base_t<num_domains>::template id<0>()].anybitset();
}

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
