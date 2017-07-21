/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/mesh/burton/burton_mesh_topology.h"
#include "flecsale/mesh/burton/burton_triangle.h"


namespace flecsale {
namespace mesh {
namespace burton {

// some type aliases
using burton_2d_triangle_t = burton_triangle_t<2>;
using burton_3d_triangle_t = burton_triangle_t<3>;

////////////////////////////////////////////////////////////////////////////////
// 2D triangle
////////////////////////////////////////////////////////////////////////////////

void burton_2d_triangle_t::update()
{
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, edge_t::domain>(this);
  const auto & a = vs[0]->coordinates();
  const auto & b = vs[1]->coordinates();
  const auto & c = vs[2]->coordinates();
  centroid_ = geom::shapes::triangle<num_dimensions>::centroid( a, b, c );
  midpoint_ = geom::shapes::triangle<num_dimensions>::midpoint( a, b, c );
  area_ = geom::shapes::triangle<num_dimensions>::area( a, b, c );
  // check the edges first
  min_length_ = abs( a - b );
  min_length_ = std::min( abs( b - c ), min_length_ );
  min_length_ = std::min( abs( c - a ), min_length_ );
}


////////////////////////////////////////////////////////////////////////////////
// 3D triangle
////////////////////////////////////////////////////////////////////////////////

void burton_3d_triangle_t::update()
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, edge_t::domain>(this);
  const auto & a = vs[0]->coordinates();
  const auto & b = vs[1]->coordinates();
  const auto & c = vs[2]->coordinates();
  centroid_ = geom::shapes::triangle<num_dimensions>::centroid( a, b, c );
  midpoint_ = geom::shapes::triangle<num_dimensions>::midpoint( a, b, c );
  area_ = geom::shapes::triangle<num_dimensions>::area( a, b, c );
  normal_ = geom::shapes::triangle<num_dimensions>::normal( a, b, c );
  // check the edges first
  min_length_ = abs( a - b );
  min_length_ = std::min( abs( b - c ), min_length_ );
  min_length_ = std::min( abs( c - a ), min_length_ );
}

} // namespace 
} // namespace 
} // namespace 
