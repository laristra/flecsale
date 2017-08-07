/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/mesh/burton/burton_wedge.h"
#include "flecsale/mesh/burton/burton_mesh_topology.h"


namespace flecsale {
namespace mesh {
namespace burton {

// some type aliases
using burton_2d_wedge_t = burton_wedge_t<2>;
using burton_3d_wedge_t = burton_wedge_t<3>;


////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 2d wedge
////////////////////////////////////////////////////////////////////////////////
template<>
void burton_2d_wedge_t::update(const mesh_topology_base_t * mesh, bool is_right)
{
  using math::abs;
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh); 
  auto vs = msh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  const auto & e = es.front()->midpoint();
  const auto & v = vs.front()->coordinates();
  if ( is_right )
    facet_normal_ = { e[1] - v[1], v[0] - e[0] };
  else
    facet_normal_ = { v[1] - e[1], e[0] - v[0] };
  facet_area_ = abs(facet_normal_);
  facet_normal_ /= facet_area_;
  facet_centroid_ = 0.5 * ( e + v );
  set_boundary( es.front()->is_boundary() );
}

////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 3d wedge
////////////////////////////////////////////////////////////////////////////////
template<>
void burton_3d_wedge_t::update(const mesh_topology_base_t * mesh, bool is_right)
{
  using math::abs;
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh); 
  auto vs = msh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  auto fs = msh->template entities<face_t::dimension, domain, face_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  assert( fs.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
  auto f = fs.front()->midpoint();
  if ( is_right )
    facet_normal_ = geom::shapes::triangle<num_dimensions>::normal( v, f, e );
  else 
    facet_normal_ = geom::shapes::triangle<num_dimensions>::normal( v, e, f );
  facet_area_ = abs(facet_normal_);
  facet_normal_ /= facet_area_;
  facet_centroid_ = geom::shapes::triangle<num_dimensions>::centroid( v, f, e );
  set_boundary( fs.front()->is_boundary() );
}

} // namespace burton
} // namespace mesh
} // namespace flecsale

