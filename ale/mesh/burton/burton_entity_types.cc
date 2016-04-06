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
#include "flecsi/mesh/mesh_topology.h"

#include "../../geom/area.h"
#include "../../geom/centroid.h"
#include "../../math/math.h"
#include "../../mesh/burton/burton_types.h"

namespace ale {
namespace mesh {

//! Type for storing instance of template specialized low level mesh.
using mesh_t = flecsi::mesh_topology_t< burton_mesh_types_t<2> >;

// aliases for ease of use
using edge_t = burton_edge_t<2>;
using cell_t = burton_cell_t<2>;
using triangle_cell_t = burton_triangle_cell_t<2>;
using quadrilateral_cell_t = burton_quadrilateral_cell_t<2>;
using polygonal_cell_t = burton_polygonal_cell_t<2>;


////////////////////////////////////////////////////////////////////////////////
// burton_edge_t
////////////////////////////////////////////////////////////////////////////////
template<>
edge_t::point_t edge_t::midpoint() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  return {0.5*(vs[0]->coordinates() + vs[1]->coordinates())};
} // burton_edge_t::midpoint

template<>
edge_t::real_t edge_t::length() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  auto & a = vs[0]->coordinates();
  auto & b = vs[1]->coordinates();
    
  using math::sqr;

  return std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
} // burton_edge_t::length

template<>
edge_t::vector_t edge_t::normal() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  auto & a = vs[0]->coordinates();
  auto & b = vs[1]->coordinates();

  return { a[1] - b[1], b[0] - a[0] };
} // burton_edge_t::normal

////////////////////////////////////////////////////////////////////////////////
// burton_cell_t
////////////////////////////////////////////////////////////////////////////////
template<>
cell_t::real_t cell_t::min_length() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  // get one of the edges as a reference
  auto edges = msh->entities<1,0>(this);
  auto min_length = edges.front()->length();
 
  for ( auto vi : vs ) {
    auto pi = vi->coordinates();
    for ( auto vj : vs ) {
      if ( vi == vj ) continue;
      auto pj = vj->coordinates();
      auto delta = pi - pj;
      min_length = std::min( abs(delta), min_length );
    }
  }
  return min_length;
}

////////////////////////////////////////////////////////////////////////////////
// burton_triangle_cell_t
////////////////////////////////////////////////////////////////////////////////

template<>
triangle_cell_t::point_t triangle_cell_t::centroid() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  return math::average( 
    vs[0]->coordinates(), vs[1]->coordinates(), vs[2]->coordinates() );
} // burton_triangle_cell_t::centroid

template<>
triangle_cell_t::real_t triangle_cell_t::area() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  auto u = vs[1]->coordinates() - vs[0]->coordinates();
  auto v = vs[2]->coordinates() - vs[0]->coordinates();
  auto cross = cross_product( u, v );
  return std::abs( cross ) / 2;
} // burton_triangle_cell_t::area
    

template<>
triangle_cell_t::real_t triangle_cell_t::min_length() const
{
  auto msh = static_cast<const mesh_t *>(mesh());

  // check the edges first
  auto edges = msh->entities<1,0>(this);
  auto eit = edges.begin();
  auto min_length = eit->length();
  std::for_each( ++eit, edges.end(), [&](auto && e) 
                 { 
                   min_length = std::min( e->length(), min_length );
                 });

  return min_length;
}

////////////////////////////////////////////////////////////////////////////////
// burton_quadrilateral_cell_t
////////////////////////////////////////////////////////////////////////////////


template<>
quadrilateral_cell_t::point_t quadrilateral_cell_t::centroid() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  return geom::centroid( vs[0]->coordinates(), vs[1]->coordinates(), 
                         vs[2]->coordinates(), vs[3]->coordinates() );
} // burton_quadrilateral_cell_t::centroid

template<>
quadrilateral_cell_t::real_t quadrilateral_cell_t::area() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  return geom::area( vs[0]->coordinates(), vs[1]->coordinates(), 
                     vs[2]->coordinates(), vs[3]->coordinates() );
} // burton_quadrilateral_cell_t::area
    

template<>
quadrilateral_cell_t::real_t quadrilateral_cell_t::min_length() const
{
  auto msh = static_cast<const mesh_t *>(mesh());

  // check the edges first
  auto edges = msh->entities<1,0>(this);
  auto eit = edges.begin();
  auto min_length = eit->length();
  std::for_each( ++eit, edges.end(), [&](auto && e) 
                 { 
                   min_length = std::min( e->length(), min_length );
                 });

  // now check the diagonal
  auto vs = msh->entities<0,0>(this);
  min_length = std::min( abs( vs[0]->coordinates() - vs[2]->coordinates() ), min_length );
  min_length = std::min( abs( vs[1]->coordinates() - vs[3]->coordinates() ), min_length );

  return min_length;
}

////////////////////////////////////////////////////////////////////////////////
// burton_polygonal_cell_t
////////////////////////////////////////////////////////////////////////////////

template<>
polygonal_cell_t::point_t polygonal_cell_t::centroid() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  std::vector<point_t> coords;
  coords.reserve( 8 );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );

  return geom::centroid( coords );
} // burton_polygonal_cell_t::centroid

template<>
polygonal_cell_t::real_t polygonal_cell_t::area() const
{
  auto msh = static_cast<const mesh_t *>(mesh());
  auto vs = msh->entities<0,0>(this);

  std::vector<point_t> coords;
  coords.reserve( 8 );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );

  return geom::area( coords );
} // burton_polygonal_cell_t::area

} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
