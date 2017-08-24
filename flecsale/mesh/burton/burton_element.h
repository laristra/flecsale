/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines a base element type.
////////////////////////////////////////////////////////////////////////////////

#pragma once


// user includes
#include "flecsale/geom/shapes/geometric_shapes.h"
#include "flecsale/geom/shapes/hexahedron.h"
#include "flecsale/geom/shapes/polygon.h"
#include "flecsale/geom/shapes/polyhedron.h"
#include "flecsale/geom/shapes/quadrilateral.h"
#include "flecsale/geom/shapes/tetrahedron.h"
#include "flecsale/geom/shapes/triangle.h"
#include "flecsale/math/general.h"
#include "flecsale/mesh/burton/burton_config.h"
#include "flecsale/utils/errors.h"
#include "flecsi/topology/mesh_types.h"

namespace flecsale {
namespace mesh {
namespace burton {

  
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

////////////////////////////////////////////////////////////////////////////////
// the list of actual coordinates
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TOPOLOGY, typename E >
auto coordinates( const MESH_TOPOLOGY * mesh, const E * ent )
{
  auto vs = mesh->template entities<0, 0>(ent);
  typename E::point_list_t coords;
  coords.reserve( vs.size() );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );
  return coords;
}
} // namespace detail



// forward decares
template< std::size_t N >
class burton_corner_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief This type provides the base for all the primal entities except 
//!        for vertices.
//!
//! \tparam NUM_DIMS  The total number of dimensions.
//! \tparam DIM       The actual dimension of the entity.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t NUM_DIMS, std::size_t DIM  >
struct burton_element_t { };

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        two-dimensional mesh edges.
//! \remark this is a 2d edge
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<2,1> :
    public flecsi::topology::mesh_entity_t<1, burton_config_t<2>::num_domains>
{

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  // the constructor
  burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the list of vertices
  auto vertices() const;

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; }

  //! \brief in 2d, this doubles as a face, so the centroid is the same as
  //!    the midpoint
  //! \remark this is only enabled in 2d
  const auto & centroid() const
  { return midpoint_; }

  //! the edge length
  auto length() const
  { return length_; }

  //! in 2d, this doubles as a face, so the area is the same as the length
  //! \remark this is only enabled in 2d
  auto area() const
  { return length_; }

  //! the edge normal
  const auto & normal() const
  { return normal_; }

  //! is this a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //! \brief update the mesh geometry
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    using math::sqr;
    using math::normal;
    auto vs = mesh->template entities<0, domain>(this);
    auto cs = mesh->template entities<2, domain>(this);
    const auto & a = vs[0]->coordinates();
    const auto & b = vs[1]->coordinates();
    midpoint_[0] = 0.5*(a[0] + b[0]);
    midpoint_[1] = 0.5*(a[1] + b[1]);
    length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
    normal_ = normal( b, a );
    normal_ /= length_;
    flags_.set( config_t::bits::boundary, (cs.size() == 1) );
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:

  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

  //! the geometric information
  real_t length_ = 0;
  point_t midpoint_ = 0;
  vector_t normal_ = 0;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh edges.
//! \remark this is a 3d edge
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,1> :
    public flecsi::topology::mesh_entity_t<1, burton_config_t<3>::num_domains>
{

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! a bitfield type
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  // the constructor
  burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; }

  //! the edge length
  auto length() const
  { return length_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! is this a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //! \brief update the mesh geometry
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    using math::sqr;
    auto vs = mesh->template entities<0, domain>(this);
    const auto & a = vs[0]->coordinates();
    const auto & b = vs[1]->coordinates();
    midpoint_[0] = 0.5*(a[0] + b[0]);
    midpoint_[1] = 0.5*(a[1] + b[1]);
    midpoint_[2] = 0.5*(a[2] + b[2]);
    length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]) );
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

  //! the geometric information
  real_t length_ = 0;
  point_t midpoint_ = 0;

}; // struct burton_edge_t

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        multi-dimensional mesh edges.
//! \tparam N The total number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_edge_t = burton_element_t<N,1>;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        two-dimensional mesh cells.
//! \remark this is a 2d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<2,2>
  : public flecsi::topology::mesh_entity_t<2, burton_config_t<2>::num_domains>
{
  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the shape type
  using shape_t = config_t::shape_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {};

  // Destructor
  ~burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  const auto & centroid() const 
  { return centroid_; };

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; };

  //! the area of the element
  auto area() const
  { return area_; };

  //! the area of the element
  auto volume() const
  { return area_; };

  //! the minimum length in the element
  auto min_length() const
  { return min_length_; }

  //! set the region id
  auto & region()
  { return region_; }

  //! get the region id
  auto region() const
  { return region_; }


  //! the element type
  auto shape() const
  { return shape_; };


  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t>
  create_entities(
    const id_t & cell,
    size_t dim,
    const connectivity_t& conn,
    id_t * entities 
  ) {
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, /* dimension */ 0 );
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

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  static 
  std::vector<size_t>
  create_bound_entities(
    size_t from_domain, 
    size_t to_domain, 
    size_t dim, 
    const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities
  ) {

    auto verts = primal_conn.get_entity_vec( cell, /* dimension */ 0 );
    auto edges = primal_conn.get_entity_vec( cell, /* dimension */ 1 );
    auto num_vertices = verts.size();

    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      // The right edge is always first
    case 0: {

      auto vp = num_vertices - 1;
      for ( auto vn=0, ind=0; vn<num_vertices; vn++ ) {
        entities[ ind++ ] = verts[vn]; // vertex 0
        entities[ ind++ ] = edges[vn]; // edge 0, abuts vertex 0
        entities[ ind++ ] = edges[vp]; // edge 3, abuts vertex 0
        vp = vn;
      }
      return std::vector<size_t>(num_vertices, 3);
    }
      //------------------------------------------------------------------------
      // wedges
      // The right wedge is always first
    case 1: {

      auto corners = domain_conn.get_entity_vec( cell, /* dimension */ 0 );

      auto vp = num_vertices - 1;
      for ( auto vn=0, ind=0; vn<num_vertices; vn++ ) {
        // wedge 0
        entities[ ind++ ] =   verts[vn]; // vertex 0
        entities[ ind++ ] =   edges[vn]; // edge 0, abuts vertex 0
        entities[ ind++ ] = corners[vn]; // corner 0
        // wedge 1
        entities[ ind++ ] =   verts[vn]; // vertex 0
        entities[ ind++ ] =   edges[vp]; // edge 3, abuts vertex 0
        entities[ ind++ ] = corners[vn]; // corner 0
        vp = vn;
      }
      return std::vector<size_t>(2*num_vertices, 3);
    }
      //------------------------------------------------------------------------
      // Failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch

  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    // get general entity connectivity
    auto vs = mesh->template entities<0, domain>(this);
    auto es = mesh->template entities<1, domain>(this);

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
        auto coords = detail::coordinates(mesh, this);
        centroid_ = geom::shapes::polygon<num_dimensions>::centroid( coords );
        midpoint_ = geom::shapes::polygon<num_dimensions>::midpoint( coords );
        area_ = geom::shapes::polygon<num_dimensions>::area( coords );
        // now check min edge length
        auto vs = mesh->template entities<0, domain>(this);
        min_length_ = detail::min_length( vs );
        break;
      }

      // there should be no unmatched case
      default:
        raise_runtime_error( "Unknown cell type" );

    } // switch
  }

private:

  //============================================================================
  // Private Data
  //============================================================================
  
  //! the region tag
  size_t region_ = 0;

  //! the geometry informtation
  real_t area_ = 0;
  point_t midpoint_ = 0;
  point_t centroid_ = 0;
  real_t min_length_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh faces.
//! \remark this is a 3d face
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,2>
  : public flecsi::topology::mesh_entity_t<2, burton_config_t<3>::num_domains>
{
  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type.
  using tag_list_t = typename config_t::tag_list_t;

  //! the shape type
  using shape_t = config_t::shape_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {};

  //! Destructor
  ~burton_element_t() = default;

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! is this a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }
  
  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //! the centroid
  const auto & centroid() const
  { return centroid_; }

  //! the midpoint used in tesselating the element
  const auto & midpoint() const 
  { return midpoint_; }

  //! the normal
  const auto & normal() const 
  { return normal_; }

  //! the area of the element
  auto area() const
  { return area_; }

  //! the minimum length in the element
  auto min_length() const
  { return min_length_; }

  //! the element type
  auto shape() const
  { return shape_; };

  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities
  )  {
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, /* dimension */ 0 );
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

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities )
  { 
    raise_implemented_error(
      "burton_3d_face_t::create_bound_entities has not been implemented"
    );
  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    using math::abs;
  
    auto vs = mesh->template entities<0, domain>(this);
    auto es = mesh->template entities<1, domain>(this);
    auto cs = mesh->template entities<3, domain>(this);
  
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
        auto coords = coordinates( mesh );
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

  //----------------------------------------------------------------------------
  //! the list of actual coordinates
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  auto coordinates(
    const MESH_TOPOLOGY * mesh,
    bool reverse = false
  ) {
    auto vs = mesh->template entities<0, domain>(this);
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

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the entity tags
  tag_list_t tags_;
  
  //! the entity flags
  bitfield_t flags_;

  //! the geometric information
  real_t area_ = 0;
  real_t min_length_ = 0;
  vector_t normal_ = 0;
  point_t centroid_ = 0;
  point_t midpoint_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;

}; // class burton_element_t


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh faces.
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_face_t = burton_element_t<N,N-1>;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh cells.
//! \remark this is a 3d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,3>
  : public flecsi::topology::mesh_entity_t<3, burton_config_t<3>::num_domains>
{
  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_element_t<3,3>;

  //! the base corner type
  using corner_t = burton_corner_t<num_dimensions>;

  //! the shape type
  using shape_t = config_t::shape_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {}

  //! Destructor
  ~burton_element_t() = default;

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  const auto & centroid() const 
  { return centroid_; };

  //! the midpoint used in tesselating the element
  const auto & midpoint() const 
  { return midpoint_; };

  //! the area of the element
  auto volume() const
  { return volume_; };

  //! the minimum length in the element
  auto min_length() const
  { return min_length_; }


  //! set the region id
  size_t & region()
  { return region_; }

  //! get the region id
  size_t region() const
  { return region_; }

  //! the element type
  auto type() const
  { return shape_; };


  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities
  )  {
    // you should only be coming in here to build edges
    assert( dim == 1 ); 
  
    // get the cell entities
    auto cell_faces = conn.get_entity_vec( cell, /* dimension */ 2 );
    // make sure the faces exist
    assert( cell_faces.size() > 0 && "no cell faces yet" );
    
    // the list of edge pairs
    std::vector< std::pair< id_t, id_t > > cell_edges;
  
    // get the edges of each face
    for ( const auto & face : cell_faces ) { 
      // get all vertices in this face
      auto face_verts = conn.get_entity_vec( face, /* dimension */ 0 );
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

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities
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
      auto cell_verts = primal_conn.get_entity_vec( cell, /* dim */ 0 );
      auto cell_edges = primal_conn.get_entity_vec( cell, /* dim */ 1 ).vec();
      auto cell_faces = primal_conn.get_entity_vec( cell, /* dim */ 2 ).vec();
  
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
        auto vert_edges = primal_conn.get_entity_vec( vert, /* dim */ 1 ).vec(); 
        auto vert_faces = primal_conn.get_entity_vec( vert, /* dim */ 2 ).vec(); 
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
      // There is an even/odd ordering so we 
      // know which way to compute normals.  So edges are defined in pairs
      //
    case 1: {
  
      // get the higher dimensional entities
      auto cell_faces   = primal_conn.get_entity_vec( cell, /* dim */ 2 );
      auto cell_corners = domain_conn.get_entity_vec( cell, /* dim */ 0 );
      
      // a counter for the number of wedges
      size_t num_wedges = 0;
  
      // loop over faces
      for ( const auto & face : cell_faces ) { 
  
        // get the vertices of the face
        auto face_verts = primal_conn.get_entity_vec( face, /* dim */ 0 );
        auto ccw_face_verts = std::vector<id_t>( face_verts.begin(), face_verts.end() );
        
        // get the edges of the face
        auto face_edges = primal_conn.get_entity_vec( face, /* dim */ 1 );
  
        // get the cells of this face
        auto face_cells = primal_conn.get_entity_vec( face, /* dim */ 3 );
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
              auto verts = primal_conn.get_entity_vec( e, /* dim */ 0 );
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
          auto point_corners = domain_conn.get_entity_vec(  *p1, /* dim */ 0 ).vec();
          auto cell_corners  = domain_conn.get_entity_vec( cell, /* dim */ 0 ).vec();
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

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    auto vs = mesh->template entities<0, domain>(this);
    auto fs = mesh->template entities<2, domain>(this);
  
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
          auto cs = mesh->template entities<3, domain>(f);
          auto reverse = (cs[0] != this); // FIXME: reverse
          auto coords = f->coordinates( mesh, reverse );
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

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the region tag
  size_t region_ = 0;

  //! the geometry informtation
  real_t volume_ = 0;
  point_t centroid_ = 0;
  point_t midpoint_ = 0;
  real_t min_length_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;

}; // class burton_element_t



////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh cells.
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_cell_t = burton_element_t<N,N>;

} // namespace burton
} // namespace mesh
} // namespace flecsale
