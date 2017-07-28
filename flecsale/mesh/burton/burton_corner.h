/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This file defines the corner entity for the burton mesh.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsale/mesh/burton/burton_vertex.h"
#include "flecsale/mesh/burton/burton_element.h"


namespace flecsale {
namespace mesh {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_corner_t type provides an interface for managing and
//!   geometry and state associated with mesh corners.
//!
//! \tparam N The domain of the corner.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_corner_t
  : public flecsi::topology::mesh_entity_t<0, burton_config_t<N>::num_domains>
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
  static constexpr auto domain = 1;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_corner_t(mesh_topology_base_t & mesh) {};

  // dissallow copying
  burton_corner_t( burton_corner_t & ) = delete;
  burton_corner_t & operator=( burton_corner_t & ) = delete;

  // dissallow moving
  burton_corner_t( burton_corner_t && ) = delete;
  burton_corner_t & operator=( burton_corner_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) { }


};

} // namespace burton
} // namespace mesh
} // namespace flecsale
