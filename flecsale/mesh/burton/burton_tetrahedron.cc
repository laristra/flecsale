/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/mesh/burton/burton_mesh_topology.h"
#include "flecsale/mesh/burton/burton_tetrahedron.h"


namespace flecsale {
namespace mesh {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
// 3D tetrahedron
////////////////////////////////////////////////////////////////////////////////

void burton_tetrahedron_t::update()
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  const auto & v0 = vs[0]->coordinates();
  const auto & v1 = vs[1]->coordinates();
  const auto & v2 = vs[2]->coordinates();
  const auto & v3 = vs[3]->coordinates();
  centroid_ = 
    geom::shapes::tetrahedron::centroid( v0, v1, v2, v3 );
  midpoint_ = 
    geom::shapes::tetrahedron::midpoint( v0, v1, v2, v3 );
  volume_ = 
    geom::shapes::tetrahedron::volume( v0, v1, v2, v3 );
  min_length_ = detail::min_length( vs );
}

} // namespace
} // namespace
} // namespace
