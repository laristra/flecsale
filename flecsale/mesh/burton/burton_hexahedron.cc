/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/mesh/burton/burton_mesh_topology.h"
#include "flecsale/mesh/burton/burton_hexahedron.h"


namespace flecsale {
namespace mesh {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
// 3D hexahedron
////////////////////////////////////////////////////////////////////////////////

void burton_hexahedron_t::update()
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  const auto & v0 = vs[0]->coordinates();
  const auto & v1 = vs[1]->coordinates();
  const auto & v2 = vs[2]->coordinates();
  const auto & v3 = vs[3]->coordinates();
  const auto & v4 = vs[4]->coordinates();
  const auto & v5 = vs[5]->coordinates();
  const auto & v6 = vs[6]->coordinates();
  const auto & v7 = vs[7]->coordinates();
  centroid_ = 
    geom::shapes::hexahedron::centroid( v0, v1, v2, v3, v4, v5, v6, v7 );
  midpoint_ = 
    geom::shapes::hexahedron::midpoint( v0, v1, v2, v3, v4, v5, v6, v7 );
  volume_ = 
    geom::shapes::hexahedron::volume( v0, v1, v2, v3, v4, v5, v6, v7 );
  min_length_ = detail::min_length( vs );
}

} // namespace
} // namespace
} // namespace
