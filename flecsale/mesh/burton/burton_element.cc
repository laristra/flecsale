/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsale/geom/shapes/hexahedron.h"
#include "flecsale/geom/shapes/polygon.h"
#include "flecsale/geom/shapes/polyhedron.h"
#include "flecsale/geom/shapes/quadrilateral.h"
#include "flecsale/geom/shapes/tetrahedron.h"
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
// the list of actual coordinates
////////////////////////////////////////////////////////////////////////////////
template<
  typename E,
  typename = std::enable_if_t< !std::is_same< E, burton_3d_face_t >::value >
>
auto coordinates(
  const typename E::mesh_topology_base_t * mesh_base,
  const E * ent
) {
  using mesh_topology_t = burton_mesh_topology_t< E::num_dimensions >;
  using vertex_t = typename E::vertex_t;
  auto mesh = static_cast<const mesh_topology_t *>(mesh_base); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(ent);
  typename E::point_list_t coords;
  coords.reserve( vs.size() );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );
  return coords;
}


////////////////////////////////////////////////////////////////////////////////
// 2d - edge
////////////////////////////////////////////////////////////////////////////////
  
// the edge midpoint
void burton_2d_edge_t::update( const mesh_topology_base_t * mesh_base )
{
  using math::sqr;
  using math::normal;
  using cell_t = burton_2d_cell_t;
  auto mesh = static_cast<const burton_2d_mesh_topology_t *>(mesh_base);
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto cs = mesh->template entities<cell_t::dimension, cell_t::domain>(this);
  const auto & a = vs[0]->coordinates();
  const auto & b = vs[1]->coordinates();
  midpoint_[0] = 0.5*(a[0] + b[0]);
  midpoint_[1] = 0.5*(a[1] + b[1]);
  length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
  normal_ = normal( b, a );
  normal_ /= length_;
  flags_.set( config_t::bits::boundary, (cs.size() == 1) );
}

////////////////////////////////////////////////////////////////////////////////
// 3d - edge
////////////////////////////////////////////////////////////////////////////////

// the edge midpoint
void burton_3d_edge_t::update( const mesh_topology_base_t * mesh_base ) 
{
  using math::sqr;
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_base);
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  const auto & a = vs[0]->coordinates();
  const auto & b = vs[1]->coordinates();
  midpoint_[0] = 0.5*(a[0] + b[0]);
  midpoint_[1] = 0.5*(a[1] + b[1]);
  midpoint_[2] = 0.5*(a[2] + b[2]);
  length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]) );
}

////////////////////////////////////////////////////////////////////////////////
// 2d - Planar Cell
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// update the cell geometry
//------------------------------------------------------------------------------
void burton_2d_cell_t::update( const mesh_topology_base_t * mesh_base )
{
  // get general entity connectivity
  auto msh = static_cast<const burton_2d_mesh_topology_t *>(mesh_base);
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
      auto coords = coordinates(msh, this);
      centroid_ = geom::shapes::polygon<num_dimensions>::centroid( coords );
      midpoint_ = geom::shapes::polygon<num_dimensions>::midpoint( coords );
      area_ = geom::shapes::polygon<num_dimensions>::area( coords );
      // now check min edge length
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
auto coordinates(
  const burton_3d_face_t::mesh_topology_base_t * mesh_base,
  const burton_3d_face_t * face,
  bool reverse = false
) {
  using vertex_t = burton_3d_face_t::vertex_t;
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_base); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(face);
  burton_3d_face_t::point_list_t coords;
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

//------------------------------------------------------------------------------
// Update goemetry
//------------------------------------------------------------------------------
void burton_3d_face_t::update( const mesh_topology_base_t * mesh_base )
{
  using math::abs;
  using cell_t = burton_3d_cell_t;

  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_base); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto es = mesh->template entities<edge_t::dimension, edge_t::domain>(this);
  auto cs = mesh->template entities<cell_t::dimension, cell_t::domain>(this);

  switch (shape_) {

    // the element is a triangle
    case shape_t::triangle: {
      const auto & a = vs[0]->coordinates();
      const auto & b = vs[1]->coordinates();
      const auto & c = vs[2]->coordinates();
      centroid_ = geom::shapes::triangle<num_dimensions>::centroid( a, b, c );
      midpoint_ = geom::shapes::triangle<num_dimensions>::midpoint( a, b, c );
      normal_ = geom::shapes::triangle<num_dimensions>::normal( a, b, c );
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
      normal_ = 
        geom::shapes::quadrilateral<num_dimensions>::normal( a, b, c, d );
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
      auto coords = coordinates( mesh_base, this );
      centroid_ = geom::shapes::polygon<num_dimensions>::centroid( coords );
      midpoint_ = geom::shapes::polygon<num_dimensions>::midpoint( coords );
      normal_ = geom::shapes::polygon<num_dimensions>::normal( coords );
      min_length_ = detail::min_length( vs );
      break;
    }

    // there should be no unmatched case
    default:
      raise_runtime_error( "Unknown cell type" );

  } // switch
  
  area_ = abs( normal_ );
  normal_ /= area_;
  flags_.set( config_t::bits::boundary, (cs.size() == 1) );
}

//------------------------------------------------------------------------------
// create_entities function for a polygon
//------------------------------------------------------------------------------
std::vector<size_t>
burton_3d_face_t::create_entities(
  const burton_3d_face_t::id_t & cell,
  size_t dim,
  const burton_3d_face_t::connectivity_t& conn,
  burton_3d_face_t::id_t * entities
)  {
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

////////////////////////////////////////////////////////////////////////////////
// 3d - Cell
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// Update goemetry
//------------------------------------------------------------------------------
void burton_3d_cell_t::update( const mesh_topology_base_t * mesh_base )
{
  auto mesh = static_cast<const burton_3d_mesh_topology_t *>(mesh_base); 
  auto vs = mesh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  auto fs = mesh->template entities<  face_t::dimension,   face_t::domain>(this);

  switch (shape_) {

    // the element is a tet
    case shape_t::tetrahedron: {
      const auto & v0 = vs[0]->coordinates();
      const auto & v1 = vs[1]->coordinates();
      const auto & v2 = vs[2]->coordinates();
      const auto & v3 = vs[3]->coordinates();
      centroid_ = 
        geom::shapes::tetrahedron::centroid( v0, v1, v2, v3 );
      midpoint_ = 
        geom::shapes::tetrahedron::midpoint( v0, v1, v2, v3 );
      volume_ = 
        geom::shapes::tetrahedron::volume( v0, v1, v2, v3 );
      min_length_ = detail::min_length( vs );
      break;
    }

    // the element is a hex
    case shape_t::hexahedron: {
      const auto & v0 = vs[0]->coordinates();
      const auto & v1 = vs[1]->coordinates();
      const auto & v2 = vs[2]->coordinates();
      const auto & v3 = vs[3]->coordinates();
      const auto & v4 = vs[4]->coordinates();
      const auto & v5 = vs[5]->coordinates();
      const auto & v6 = vs[6]->coordinates();
      const auto & v7 = vs[7]->coordinates();
      centroid_ = 
        geom::shapes::hexahedron::centroid( v0, v1, v2, v3, v4, v5, v6, v7 );
      midpoint_ = 
        geom::shapes::hexahedron::midpoint( v0, v1, v2, v3, v4, v5, v6, v7 );
      volume_ = 
        geom::shapes::hexahedron::volume( v0, v1, v2, v3, v4, v5, v6, v7 );
      min_length_ = detail::min_length( vs );
      break;
    }

    // the element is a polyhedron
    case shape_t::polyhedron: {
      geom::shapes::polyhedron<point_t> poly;     
      for ( auto f : fs ) {
        auto cs = mesh->template entities<cell_t::dimension, cell_t::domain>(f);
        auto reverse = (cs[0] != this); // FIXME: reverse
        auto coords = coordinates( mesh_base, f, reverse );
        poly.insert( coords );
      }
      centroid_ = poly.centroid();
      midpoint_ = poly.midpoint();
      volume_ = poly.volume();
      min_length_ = detail::min_length( vs );
      break;
    }

    // there should be no unmatched case
    default:
      raise_runtime_error( "Unknown cell type" );

  } // switch
  
}

//------------------------------------------------------------------------------
// create_entities function for a polyhedron
//------------------------------------------------------------------------------
std::vector<size_t>
burton_3d_cell_t::create_entities(
  const burton_3d_cell_t::id_t & cell,
  size_t dim,
  const burton_3d_cell_t::connectivity_t& conn,
  burton_3d_cell_t::id_t * entities
)  {
  // you should only be coming in here to build edges
  assert( dim == 1 ); 

  // get the cell entities
  auto cell_faces = conn.get_entity_vec( cell, face_t::dimension );
  // make sure the faces exist
  assert( cell_faces.size() > 0 && "no cell faces yet" );
  
  // the list of edge pairs
  std::vector< std::pair< id_t, id_t > > cell_edges;

  // get the edges of each face
  for ( const auto & face : cell_faces ) { 
    // get all vertices in this face
    auto face_verts = conn.get_entity_vec( face, vertex_t::dimension );
    assert( face_verts.size() > 2 && "not enough vertices for a valid face" );
    // reserve space 
    cell_edges.reserve( cell_edges.size() + face_verts.size() );
    // add each edge pair to the list
    auto vp = std::prev( face_verts.end() );
    assert( vp != face_verts.begin() && "no vertices in this face" ); 
    for ( auto vn = face_verts.begin(); vn != face_verts.end(); vn++ ) {
      assert( *vn != *vp && "edge has two equal vertices" );
      if ( *vn < *vp ) 
        cell_edges.emplace_back( std::make_pair( *vn, *vp ) );
      else
        cell_edges.emplace_back( std::make_pair( *vp, *vn ) );
      vp = vn;
    }          
  }

  // now sort the list of edges
  std::sort( 
    cell_edges.begin(), cell_edges.end(), 
    [](const auto & a, const auto & b) 
    { 
      id_t a1 = a.first;
      id_t b1 = b.first;
      if ( a1 == b1 )
        return ( a.second < b.second );
      else
        return (a1 < b1);
    }
  );
  // remove uniques
  auto end = std::unique( cell_edges.begin(), cell_edges.end() );
  auto num_edges = std::distance( cell_edges.begin(), end );

  // copy the unique results to the output array   
  size_t i = 0;

  std::for_each( 
    cell_edges.begin(), end, 
    [&](const auto & edge) 
    {
      entities[i++] = edge.first;
      entities[i++] = edge.second;
    }
  );

  return std::vector<size_t>( num_edges, 2 );
}

//------------------------------------------------------------------------------
// create_bound_entities function for polyhedron
//------------------------------------------------------------------------------
std::vector<size_t>
burton_3d_cell_t::create_bound_entities(
  size_t from_domain,
  size_t to_domain,
  size_t dim,
  const burton_3d_cell_t::id_t & cell,
  const burton_3d_cell_t::connectivity_t& primal_conn,
  const burton_3d_cell_t::connectivity_t& domain_conn,
  burton_3d_cell_t::id_t * entities
) {

  size_t i = 0;

  switch (dim) {
    //------------------------------------------------------------------------
    // Corners
    //
    // Take your right hand, its origin is the vertex of the corner.  Curl 
    // your hand from the first edge to the second edge, with the third edge
    // aligned with your thumb.  You hand also curls from the first to the 
    // first to second face, with the third face on the bottom.
    //
  case 0: {

    std::vector<size_t> entity_count;

    // get the cell entities
    auto cell_verts = primal_conn.get_entity_vec( cell, vertex_t::dimension );
    auto cell_edges = primal_conn.get_entity_vec( cell, edge_t::dimension ).vec();
    auto cell_faces = primal_conn.get_entity_vec( cell, face_t::dimension ).vec();

    // sort the edges and faces for intersections later      
    std::sort( cell_edges.begin(), cell_edges.end() );
    std::sort( cell_faces.begin(), cell_faces.end() );

    // temparary lists
    std::vector<id_t> edges, faces;

    for ( const auto & vert : cell_verts ) { 
      // clear temporary lits
      edges.clear();
      faces.clear();
      // get all entities attached to this vertex
      auto vert_edges = primal_conn.get_entity_vec( vert, edge_t::dimension ).vec(); 
      auto vert_faces = primal_conn.get_entity_vec( vert, face_t::dimension ).vec(); 
      // sort the lists for intersectinos
      std::sort( vert_edges.begin(), vert_edges.end() );
      std::sort( vert_faces.begin(), vert_faces.end() );
      // get the intersections of the sets
      std::set_intersection( vert_edges.begin(), vert_edges.end(),
                             cell_edges.begin(), cell_edges.end(),
                             std::back_inserter(edges));
      std::set_intersection( vert_faces.begin(), vert_faces.end(),
                             cell_faces.begin(), cell_faces.end(),
                             std::back_inserter(faces));
      // add the entities to the list
      entities[i++] = vert;
      for ( auto e : edges ) entities[i++] = e;
      for ( auto f : faces ) entities[i++] = f;
      // add the final number of entities to the list
      entity_count.emplace_back( 1 + edges.size() + faces.size() );
      
    }  // foreach vertex

    return entity_count;
  }
    //------------------------------------------------------------------------
    // Wedges
    //
    // Each corner has 6 vertices.  There is an even/odd ordering so we 
    // no which way to compute normals.  So edges are defined in pairs
    //
  case 1: {

    // get the higher dimensional entities
    auto cell_faces   = primal_conn.get_entity_vec( cell,   face_t::dimension );
    auto cell_corners = domain_conn.get_entity_vec( cell, corner_t::dimension );
    
    // a counter for the number of wedges
    size_t num_wedges = 0;

    // loop over faces
    for ( const auto & face : cell_faces ) { 

      // get the vertices of the face
      auto face_verts = primal_conn.get_entity_vec( face, vertex_t::dimension );
      auto ccw_face_verts = std::vector<id_t>( face_verts.begin(), face_verts.end() );
      
      // get the edges of the face
      auto face_edges = primal_conn.get_entity_vec( face, edge_t::dimension );

      // get the cells of this face
      auto face_cells = primal_conn.get_entity_vec( face, cell_t::dimension );
      // reverse the list of vertices if this is backwards
      if ( face_cells[0] != cell ) 
        std::reverse( ccw_face_verts.begin(), ccw_face_verts.end() );

      // a lambda function to locate an edge connected to two points
      auto _find_edge = [&]( const auto & pa, const auto & pb  ) 
      {
        // locate the edge with the two vertices
        auto edge = std::find_if( 
          face_edges.begin(), face_edges.end(), 
          [&]( const auto & e ) 
          { 
            auto verts = primal_conn.get_entity_vec( e, vertex_t::dimension );
            assert( verts.size() == 2 && "should be two vertices per edge" );
            return ( (verts[0] == pa && verts[1] == pb) || 
                     (verts[0] == pb && verts[1] == pa) );
          } 
        );
        // make sure we found an edge
        assert( edge != face_edges.end() );
        // return the edge
        return edge;
      };
        
      // loop over each edge (pair of vertices of the face)
      // there is a wedge for each vertex -> edge -> face combination
      auto p1 = std::prev( ccw_face_verts.end() );
      auto p0 = std::prev( p1 );
      auto edge0 = _find_edge( *p0, *p1 );
      for ( auto p2=ccw_face_verts.begin(); p2!=ccw_face_verts.end(); p2++ ) {
        // get the next edge
        auto edge1 = _find_edge( *p1, *p2 );
        // get the corner of this point, but also associated with this cell
        auto point_corners = domain_conn.get_entity_vec(  *p1, corner_t::dimension ).vec();
        auto cell_corners  = domain_conn.get_entity_vec( cell, corner_t::dimension ).vec();
        // sort the lists for intersectinos
        std::sort( point_corners.begin(), point_corners.end() );
        std::sort( cell_corners.begin(),  cell_corners.end() );
        // get the intersections of the sets
        std::vector<id_t> corners;
        corners.reserve(1);
        std::set_intersection( point_corners.begin(), point_corners.end(),
                               cell_corners.begin(), cell_corners.end(),
                               std::back_inserter(corners));
        assert( corners.size() == 1 );
        const auto & c1 = corners.front();
        // the first point, this is the right (even) one
        entities[i++] = *p1;
        entities[i++] = *edge0;
        entities[i++] = face;
        entities[i++] = c1;
        // for the next point, this one is the left (odd) one
        entities[i++] = *p1;
        entities[i++] = *edge1;
        entities[i++] = face;
        entities[i++] = c1;
        // update the iterators
        p0 = p1;
        p1 = p2;
        edge0 = edge1;
        // increment the wedge counter
        num_wedges += 2;
      }
    } // foreach vertex

    return std::vector<size_t>( num_wedges, 4 );
  }
    //------------------------------------------------------------------------
    // failure
  default:
    raise_runtime_error("Unknown bound entity type");
    return {};

  } // switch
}

} // namespace
} // namespace
} // namespace
