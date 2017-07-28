/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/geom/shapes/polygon.h"
#include "flecsale/geom/shapes/quadrilateral.h"
#include "flecsale/geom/shapes/triangle.h"
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


namespace detail {

////////////////////////////////////////////////////////////////////////////////
/// \brief compute the minimum length
////////////////////////////////////////////////////////////////////////////////
template< typename VS >
inline auto min_length( VS && vs ) {

  std::decay_t< decltype(vs[0]->coordinates()[0]) > len;
  bool first = true;

  // check each vertex combination
  for ( auto vi : vs ) {
    const auto & pi = vi->coordinates();
    for ( auto vj : vs ) {
      if ( vi == vj ) continue;
      const auto & pj = vj->coordinates();
      auto delta = pi - pj;
      if ( first ) {
        len = abs(delta);
        first = false;
      }
      else {
        len = std::min( abs(delta), len );
      }
    }
  }

  return len;
}
} // namespace detail


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

//------------------------------------------------------------------------------
// update the cell geometry
//------------------------------------------------------------------------------
void burton_2d_cell_t::update()
{
  // get general entity connectivity
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto es = msh->template entities<edge_t::dimension, edge_t::domain>(this);

  switch (shape_) {

    // the element is a triangle
    case shape_t::triangle: {
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
      break;
    }

    // the element is a quadrilateral
    case shape_t::quadrilateral: {
      const auto & a = vs[0]->coordinates();
      const auto & b = vs[1]->coordinates();
      const auto & c = vs[2]->coordinates();
      const auto & d = vs[3]->coordinates();
      centroid_ = 
        geom::shapes::quadrilateral<num_dimensions>::centroid( a, b, c, d );
      midpoint_ = 
        geom::shapes::quadrilateral<num_dimensions>::midpoint( a, b, c, d );
      area_ = 
        geom::shapes::quadrilateral<num_dimensions>::area( a, b, c, d );
      // check the edges first
      min_length_ = abs( a - b );
      min_length_ = std::min( abs( b - c ), min_length_ );
      min_length_ = std::min( abs( c - d ), min_length_ );
      min_length_ = std::min( abs( d - a ), min_length_ );
      // now check the diagonal
      min_length_ = std::min( abs( a - c ), min_length_ );
      min_length_ = std::min( abs( b - d ), min_length_ );
      break;
    }

    // the element is a polygon
    case shape_t::polygon: {
      auto coords = coordinates();
      centroid_ = geom::shapes::polygon<num_dimensions>::centroid( coords );
      midpoint_ = geom::shapes::polygon<num_dimensions>::midpoint( coords );
      area_ = geom::shapes::polygon<num_dimensions>::area( coords );
      // now check min edge length
      auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh()); 
      auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
      min_length_ = detail::min_length( vs );
      break;
    }

    // there should be no unmatched case
    default:
      raise_runtime_error( "Unknown cell type" );

  } // switch

}

//------------------------------------------------------------------------------
// create_entities function for polygon
//------------------------------------------------------------------------------
std::vector<size_t> 
burton_2d_cell_t::create_entities(
  const id_t & cell, 
  size_t dim,
  const burton_2d_cell_t::connectivity_t& conn,
  burton_2d_cell_t::id_t * entities
) {
  
  assert( dim == 1 );

  auto v = conn.get_entity_vec( cell, vertex_t::dimension );
  auto num_cell_verts = v.size();

  size_t ind=0;
  for ( auto i=0; i<num_cell_verts-1; i++ ) {
    auto vp = v[i];
    auto vn = v[i+1];
    entities[ ind++ ] = vp;
    entities[ ind++ ] = vn;
  }
  entities[ ind++ ] = v[ num_cell_verts-1 ];
  entities[ ind++ ] = v[ 0 ];

  return std::vector<size_t>(num_cell_verts, 2);
}

//------------------------------------------------------------------------------
// create_bound_entities function for polygonal
//------------------------------------------------------------------------------
std::vector<size_t>
burton_2d_cell_t::create_bound_entities(
  size_t from_domain,
  size_t to_domain,
  size_t dim,
  const burton_2d_cell_t::id_t & cell,
  const burton_2d_cell_t::connectivity_t& primal_conn,
  const burton_2d_cell_t::connectivity_t& domain_conn,
  burton_2d_cell_t::id_t * c
) {

  auto verts = primal_conn.get_entity_vec( cell, vertex_t::dimension );
  auto edges = primal_conn.get_entity_vec( cell,   edge_t::dimension );
  auto num_vertices = verts.size();

  switch (dim) {
    //------------------------------------------------------------------------
    // Corners
    // The right edge is always first
  case 0: {

    auto vp = num_vertices - 1;
    for ( auto vn=0, ind=0; vn<num_vertices; vn++ ) {
      c[ ind++ ] = verts[vn]; // vertex 0
      c[ ind++ ] = edges[vn]; // edge 0, abuts vertex 0
      c[ ind++ ] = edges[vp]; // edge 3, abuts vertex 0
      vp = vn;
    }
    return std::vector<size_t>(num_vertices, 3);
  }
    //------------------------------------------------------------------------
    // wedges
    // The right wedge is always first
  case 1: {

    auto corners = domain_conn.get_entity_vec( cell, corner_t::dimension );

    auto vp = num_vertices - 1;
    for ( auto vn=0, ind=0; vn<num_vertices; vn++ ) {
      // wedge 0
      c[ ind++ ] =   verts[vn]; // vertex 0
      c[ ind++ ] =   edges[vn]; // edge 0, abuts vertex 0
      c[ ind++ ] = corners[vn]; // corner 0
      // wedge 1
      c[ ind++ ] =   verts[vn]; // vertex 0
      c[ ind++ ] =   edges[vp]; // edge 3, abuts vertex 0
      c[ ind++ ] = corners[vn]; // corner 0
      vp = vn;
    }
    return std::vector<size_t>(2*num_vertices, 3);
  }
    //------------------------------------------------------------------------
    // Failure
  default:
    raise_runtime_error("Unknown bound entity type");
  } // switch

} // create_bound_entities


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
