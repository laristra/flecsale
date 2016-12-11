/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "ale/mesh/burton/burton_vertex.h"
#include "ale/mesh/burton/burton_mesh_topology.h"


namespace ale {
namespace mesh {

// some type aliases
using burton_2d_vertex_t = burton_vertex_t<2>;
using burton_3d_vertex_t = burton_vertex_t<3>;


////////////////////////////////////////////////////////////////////////////////
// is this a boundary
////////////////////////////////////////////////////////////////////////////////
bool burton_2d_vertex_t::is_boundary() const
{
  auto flags = get_accessor(*mesh_, mesh, node_flags, bitfield_t, dense, 0);
  return flags[mesh_entity_base_t<num_domains>::template id<0>()].bitset( config_t::bits::boundary );
}

bool burton_3d_vertex_t::is_boundary() const
{
  auto flags = get_accessor(*mesh_, mesh, node_flags, bitfield_t, dense, 0);
  return flags[mesh_entity_base_t<num_domains>::template id<0>()].bitset( config_t::bits::boundary );
}

////////////////////////////////////////////////////////////////////////////////
// tag as a boundary
////////////////////////////////////////////////////////////////////////////////
void burton_2d_vertex_t::tag(const burton_2d_vertex_t::tag_t & tag)
{
  auto flags = get_accessor(*mesh_, mesh, node_tags, tag_list_t, dense, 0);
  flags[mesh_entity_base_t<num_domains>::template id<0>()].push_back( tag );
}

void burton_3d_vertex_t::tag(const burton_3d_vertex_t::tag_t & tag)
{
  auto flags = get_accessor(*mesh_, mesh, node_tags, tag_list_t, dense, 0);
  flags[mesh_entity_base_t<num_domains>::template id<0>()].push_back( tag );
}

////////////////////////////////////////////////////////////////////////////////
// get the boundary tags
////////////////////////////////////////////////////////////////////////////////
const burton_2d_vertex_t::tag_list_t &  burton_2d_vertex_t::tags() const
{
  auto flags = get_accessor(*mesh_, mesh, node_tags, tag_list_t, dense, 0);
  return flags[mesh_entity_base_t<num_domains>::template id<0>()];
}

const burton_3d_vertex_t::tag_list_t &  burton_3d_vertex_t::tags() const
{
  auto flags = get_accessor(*mesh_, mesh, node_tags, tag_list_t, dense, 0);
  return flags[mesh_entity_base_t<num_domains>::template id<0>()];
}

} // namespace mesh
} // namespace ale

