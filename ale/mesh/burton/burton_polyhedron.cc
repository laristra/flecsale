/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "ale/mesh/burton/burton_mesh_topology.h"
#include "ale/mesh/burton/burton_polyhedron.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
// 3D polyhedron
////////////////////////////////////////////////////////////////////////////////

// the centroid
burton_polyhedron_t::point_t burton_polyhedron_t::centroid() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto fs = msh->template entities<  face_t::dimension,   face_t::domain>(this);
  geom::shapes::polyhedron<point_t> poly;     
  for ( auto f : fs ) {
    auto cs = msh->template entities<cell_t::dimension, cell_t::domain>(f);
    auto reverse = (cs[0] != this); // FIXME: reverse
    poly.insert( f->coordinates(reverse) );
  }
  return poly.centroid();
}

// the midpoint
burton_polyhedron_t::point_t burton_polyhedron_t::midpoint() const
{
  auto coords = coordinates();
  return geom::shapes::polyhedron<point_t>::midpoint( coords );
}


// the area of the cell
burton_polyhedron_t::real_t burton_polyhedron_t::volume() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto fs = msh->template entities<  face_t::dimension,   face_t::domain>(this);
  geom::shapes::polyhedron<point_t> poly;     
  for ( auto f : fs ) 
    poly.insert( f->coordinates() );
  return poly.volume();
}


} // namespace
} // namespace
