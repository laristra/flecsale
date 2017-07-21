/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/mesh/burton/burton_mesh_topology.h"
#include "flecsale/mesh/burton/burton_polyhedron.h"


namespace flecsale {
namespace mesh {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
// 3D polyhedron
////////////////////////////////////////////////////////////////////////////////

void burton_polyhedron_t::update()
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto fs = msh->template entities<  face_t::dimension,   face_t::domain>(this);
  geom::shapes::polyhedron<point_t> poly;     
  for ( auto f : fs ) {
    auto cs = msh->template entities<cell_t::dimension, cell_t::domain>(f);
    auto reverse = (cs[0] != this); // FIXME: reverse
    poly.insert( f->coordinates(reverse) );
  }
  centroid_ = poly.centroid();
  midpoint_ = poly.midpoint();
  volume_ = poly.volume();
}


} // namespace
} // namespace
} // namespace
