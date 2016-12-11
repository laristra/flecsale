/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "ale/mesh/burton/burton_mesh_topology.h"
#include "ale/mesh/burton/burton_quadrilateral.h"


namespace ale {
namespace mesh {

// some type aliases
using burton_2d_quad_t = burton_quadrilateral_t<2>;
using burton_3d_quad_t = burton_quadrilateral_t<3>;

////////////////////////////////////////////////////////////////////////////////
// 2D triangle
////////////////////////////////////////////////////////////////////////////////

//! the centroid
burton_2d_quad_t::point_t burton_2d_quad_t::centroid() const
{
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::centroid( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}


//! the midpoint
burton_2d_quad_t::point_t burton_2d_quad_t::midpoint() const
{
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::midpoint( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

//! the area of the cell
burton_2d_quad_t::real_t burton_2d_quad_t::area() const
{
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::area( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

//! the minimum length in the cell
burton_2d_quad_t::real_t burton_2d_quad_t::min_length() const
{
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, edge_t::domain>(this);
  // check the edges first
  auto eit = es.begin();
  auto min_length = (*eit)->length();
  std::for_each( ++eit, es.end(), [&](auto && e) 
                 { 
                   min_length = std::min( (*e)->length(), min_length );
                 });
  // now check the diagonal
  min_length = std::min( abs( vs[0]->coordinates() - vs[2]->coordinates() ), min_length );
  min_length = std::min( abs( vs[1]->coordinates() - vs[3]->coordinates() ), min_length );
  // return the result
  return min_length;
}


////////////////////////////////////////////////////////////////////////////////
// 3D triangle
////////////////////////////////////////////////////////////////////////////////

//! the centroid
burton_3d_quad_t::point_t burton_3d_quad_t::centroid() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::centroid( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

//! the midpoint
burton_3d_quad_t::point_t burton_3d_quad_t::midpoint() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::midpoint( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}


//! the normal
burton_3d_quad_t::vector_t burton_3d_quad_t::normal() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::normal( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

//! the area of the cell
burton_3d_quad_t::real_t burton_3d_quad_t::area() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geom::shapes::quadrilateral<num_dimensions>::area( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

//! the minimum length in the cell
burton_3d_quad_t::real_t burton_3d_quad_t::min_length() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, edge_t::domain>(this);
  // check the edges first
  auto eit = es.begin();
  auto min_length = (*eit)->length();
  std::for_each( ++eit, es.end(), [&](auto && e) 
                 { 
                   min_length = std::min( (*e)->length(), min_length );
                 });
  // now check the diagonal
  min_length = std::min( abs( vs[0]->coordinates() - vs[2]->coordinates() ), min_length );
  min_length = std::min( abs( vs[1]->coordinates() - vs[3]->coordinates() ), min_length );
  // return the result
  return min_length;
}


} // namespace
} // namespace
