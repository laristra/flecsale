/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "ale/mesh/burton/burton_mesh_topology.h"
#include "ale/mesh/burton/burton_tetrahedron.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
// 3D tetrahedron
////////////////////////////////////////////////////////////////////////////////

// the centroid
burton_tetrahedron_t::point_t burton_tetrahedron_t::centroid() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::tetrahedron::centroid( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

// the midpoint
burton_tetrahedron_t::point_t burton_tetrahedron_t::midpoint() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::tetrahedron::midpoint( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}


// the area of the cell
burton_tetrahedron_t::real_t burton_tetrahedron_t::volume() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::tetrahedron::volume( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

} // namespace
} // namespace
