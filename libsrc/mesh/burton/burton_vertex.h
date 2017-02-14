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
#include "geom/shapes/geometric_shapes.h"
#include "mesh/burton/burton_config.h"
#include "utils/errors.h"
#include "flecsi/topology/mesh_types.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//!
//! \tparam N The dimension of the vertex.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_vertex_t {};

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//! \remark This is a two dimensional specialization.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_vertex_t<2> : public 
  flecsi::topology::mesh_entity_t<0, burton_config_t<2>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::topology::mesh_topology_base_t;
 
  //! the mesh traits
  using config_t = burton_config_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_vertex_t(mesh_topology_base_t & mesh) : mesh_(&mesh) 
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

  //! return true if this is on a boundary
  bool is_boundary() const;

  //! get all entity tags
  const tag_list_t & tags() const;
  //! tag entity
  void tag(const tag_t & tag);
  //! does entity have a tag
  bool has_tag(const tag_t & tag);


  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  mesh_topology_base_t * mesh_ = nullptr;

  //! the coordinates of the vertex
  point_t coordinates_;


};


////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//! \remark This is a three dimensional specialization.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_vertex_t<3> : 
    public flecsi::topology::mesh_entity_t<0, burton_config_t<2>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::topology::mesh_topology_base_t;
 
  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

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
  burton_vertex_t(mesh_topology_base_t & mesh) : mesh_(&mesh) 
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

  //! return true if this is on a boundary
  bool is_boundary() const;

  //! get all entity tags
  const tag_list_t & tags() const;
  //! tag entity
  void tag(const tag_t & tag);
  //! does entity have a tag
  bool has_tag(const tag_t & tag);

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  mesh_topology_base_t * mesh_ = nullptr;
  
  //! the coordinates of the vertex
  point_t coordinates_;

}; // class burton_vertex_t



} // namespace mesh
} // namespace ale

