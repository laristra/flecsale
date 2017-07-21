/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/math/general.h"
#include "flecsale/mesh/burton/burton_mesh_topology.h"
#include "flecsale/mesh/burton/burton_element.h"


namespace flecsale {
namespace mesh {
namespace burton {

// some type aliases
using burton_2d_edge_t = burton_element_t<2,1>;
using burton_2d_cell_t = burton_element_t<2,2>;
using burton_3d_edge_t = burton_element_t<3,1>;
using burton_3d_face_t = burton_element_t<3,2>;
using burton_3d_cell_t = burton_element_t<3,3>;

using flecsi::topology::mesh_entity_base_t;


////////////////////////////////////////////////////////////////////////////////
// 2d - edge
////////////////////////////////////////////////////////////////////////////////

// the list of actual coordinates
burton_2d_edge_t::point_list_t burton_2d_edge_t::coordinates() const
{
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  point_list_t coords;
  coords.front() = vs.front()->coordinates();
  coords.back () = vs.back ()->coordinates();
  return coords;
}

  
// the edge midpoint
void burton_2d_edge_t::update()
{
  using math::sqr;
  using math::normal;
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  const auto & a = vs[0]->coordinates();
  const auto & b = vs[1]->coordinates();
  midpoint_[0] = 0.5*(a[0] + a[0]);
  midpoint_[1] = 0.5*(a[1] + b[1]);
  length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
  normal_ = normal( b, a );

}

//! is this a boundary
bool burton_2d_edge_t::is_boundary() const
{
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto cs = mesh->template entities<burton_2d_cell_t::dimension, burton_2d_cell_t::domain>(this);
  return (cs.size() == 1);
}

////////////////////////////////////////////////////////////////////////////////
// 3d - edge
////////////////////////////////////////////////////////////////////////////////

// the list of actual coordinates
burton_3d_edge_t::point_list_t burton_3d_edge_t::coordinates() const
{
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<0,0>(this); 
  point_list_t coords;
  coords.front() = vs.front()->coordinates();
  coords.back () = vs.back ()->coordinates();
  return coords;
}

  
// the edge midpoint
void burton_3d_edge_t::update() 
{
  using math::sqr;
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  const auto & a = vs[0]->coordinates();
  const auto & b = vs[1]->coordinates();
  midpoint_[0] = 0.5*(a[0] + a[0]);
  midpoint_[1] = 0.5*(a[1] + b[1]);
  midpoint_[2] = 0.5*(a[2] + b[2]);
  length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]) );
}

////////////////////////////////////////////////////////////////////////////////
// 2d - Planar Cell
////////////////////////////////////////////////////////////////////////////////


// the list of actual coordinates
burton_2d_cell_t::point_list_t burton_2d_cell_t::coordinates() const
{
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  point_list_t coords;
  coords.reserve( vs.size() );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );
  return coords;
}


////////////////////////////////////////////////////////////////////////////////
// 3d - Face
////////////////////////////////////////////////////////////////////////////////

// the list of actual coordinates
burton_3d_face_t::point_list_t burton_3d_face_t::coordinates( bool reverse ) const
{
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  point_list_t coords;
  if ( reverse ) {
    coords.resize( vs.size() );
    size_t cnt = vs.size()-1;
    for ( auto v : vs ) 
      coords[cnt--] = v->coordinates();
  }
  else {
    coords.reserve( vs.size() );
    for ( auto v : vs ) 
      coords.emplace_back( v->coordinates() );     
  }
  return coords;
}

//! is this a boundary
bool burton_3d_face_t::is_boundary() const
{
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto cs = mesh->template entities<burton_3d_cell_t::dimension, burton_3d_cell_t::domain>(this);
  return (cs.size() == 1);
}

////////////////////////////////////////////////////////////////////////////////
// 3d - Cell
////////////////////////////////////////////////////////////////////////////////

// the list of actual coordinates
burton_3d_cell_t::point_list_t burton_3d_cell_t::coordinates() const
{
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  point_list_t coords;
  coords.reserve( vs.size() );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );
  return coords;
}

} // namespace
} // namespace
} // namespace
