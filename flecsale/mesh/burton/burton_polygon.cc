/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/mesh/burton/burton_mesh_topology.h"
#include "flecsale/mesh/burton/burton_polygon.h"


namespace flecsale {
namespace mesh {
namespace burton {

// some type aliases
using burton_2d_polygon_t = burton_polygon_t<2>;
using burton_3d_polygon_t = burton_polygon_t<3>;

////////////////////////////////////////////////////////////////////////////////
// 2D polygon
////////////////////////////////////////////////////////////////////////////////


void burton_2d_polygon_t::update()
{
  auto coords = coordinates();
  centroid_ = geom::shapes::polygon<num_dimensions>::centroid( coords );
  midpoint_ = geom::shapes::polygon<num_dimensions>::midpoint( coords );
  area_ = geom::shapes::polygon<num_dimensions>::area( coords );
  // now check min edge length
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  min_length_ = detail::min_length( vs );

}


////////////////////////////////////////////////////////////////////////////////
// 3D polygon
////////////////////////////////////////////////////////////////////////////////


void burton_3d_polygon_t::update()
{
  auto coords = coordinates();
  centroid_ = geom::shapes::polygon<num_dimensions>::centroid( coords );
  midpoint_ = geom::shapes::polygon<num_dimensions>::midpoint( coords );
  area_ = geom::shapes::polygon<num_dimensions>::area( coords );
  normal_ = geom::shapes::polygon<num_dimensions>::normal( coords );
  // now check min edge length
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  min_length_ = detail::min_length( vs );
}

} // namespace
} // namespace
} // namespace
