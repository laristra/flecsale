/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Provides the vertex type for burton_mesh_t. 
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

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//!
//! \tparam N The dimension of the vertex.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_vertex_t : public 
  flecsi::topology::mesh_entity_t<0, burton_config_t<N>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<N>;

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

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type.
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  template< typename...ARGS >
  burton_vertex_t(ARGS &&... args) 
    : coordinates_{ std::forward<ARGS>(args)... }
  {}

  // dissallow copying
  burton_vertex_t( burton_vertex_t & ) = delete;
  burton_vertex_t & operator=( burton_vertex_t & ) = delete;

  // dissallow moving
  burton_vertex_t( burton_vertex_t && ) = delete;
  burton_vertex_t & operator=( burton_vertex_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  const point_t & coordinates() const noexcept
  { return coordinates_; }

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  //! \remark this is the non const version
  point_t & coordinates() noexcept
  { return coordinates_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! return true if this is on a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the coordinates of the vertex
  point_t coordinates_ = 0;

  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

}; // class burton_vertex_t



} // namespace burton
} // namespace mesh
} // namespace flecsale

