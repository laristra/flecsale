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
#include "flecsale/mesh/burton/burton_config.h"
#include "flecsale/utils/errors.h"
#include "flecsi/topology/mesh_storage.h"
#include "flecsi/topology/mesh_types.h"

namespace flecsale {
namespace mesh {
namespace burton {

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

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

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

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

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
  void update( const mesh_topology_base_t * mesh );

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

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

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

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

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
  void update( const mesh_topology_base_t * mesh ); 

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

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

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
  using id_t = flecsi::utils::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = flecsi::topology::domain_connectivity<num_dimensions>;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base corner type
  using corner_t = burton_corner_t<num_dimensions>;

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
  ); 

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
    const id_t & cell_id,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities ); 

  //! \brief update the mesh geometry
  void update( const mesh_topology_base_t * mesh ); 

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

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

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
  using id_t = flecsi::utils::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = flecsi::topology::domain_connectivity<num_dimensions>;

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
    id_t * entities ); 

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

  //! \brief update the mesh geometry
  void update( const mesh_topology_base_t * mesh ); 

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

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  // the flecsi id type
  using id_t = flecsi::utils::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = flecsi::topology::domain_connectivity<num_dimensions>;

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
    id_t * entities );

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
    id_t * entities );

  //! \brief update the mesh geometry
  void update( const mesh_topology_base_t * mesh ); 

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
