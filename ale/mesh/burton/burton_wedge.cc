/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

//! user includes
#include "ale/mesh/burton/burton_wedge.h"
#include "ale/mesh/burton/burton_mesh_topology.h"


namespace ale {
namespace mesh {

// some type aliases
using burton_2d_wedge_t = burton_wedge_t<2>;
using burton_3d_wedge_t = burton_wedge_t<3>;


////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 2d wedge
////////////////////////////////////////////////////////////////////////////////
burton_2d_wedge_t::vector_t burton_2d_wedge_t::facet_normal_left() const
{
  using flecsi::mesh_entity_base_t;
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto vs = msh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
  return { v[1] - e[1], e[0] - v[0] };
}

burton_2d_wedge_t::vector_t burton_2d_wedge_t::facet_normal_right() const
{
  using flecsi::mesh_entity_base_t;
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh_); 
  auto vs = msh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
    return { e[1] - v[1], v[0] - e[0] };
}


////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 3d wedge
////////////////////////////////////////////////////////////////////////////////
burton_3d_wedge_t::vector_t burton_3d_wedge_t::facet_normal_left() const
{
  using flecsi::mesh_entity_base_t;
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto vs = msh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  auto fs = msh->template entities<face_t::dimension, domain, face_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  assert( fs.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
  auto f = fs.front()->centroid();
  return geom::triangle<num_dimensions>::normal( v, e, f );
}

burton_3d_wedge_t::vector_t burton_3d_wedge_t::facet_normal_right() const
{
  using flecsi::mesh_entity_base_t;
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh_); 
  auto vs = msh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  auto fs = msh->template entities<face_t::dimension, domain, face_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  assert( fs.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
  auto f = fs.front()->centroid();
  return geom::triangle<num_dimensions>::normal( v, f, e );
}

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
