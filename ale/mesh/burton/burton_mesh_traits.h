/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines the mesh traits for the burton_mesh_t class.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "ale/common/types.h"
#include "ale/geom/point.h"
#include "ale/math/vector.h"
#include "ale/utils/fixed_vector.h"
#include "flecsi/utils/bitfield.h"
#include "flecsi/data/data.h"

namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
/// \brief The mesh traits
/// \tparam N  The number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_mesh_traits_t {

  //! the size type
  using size_t = std::size_t;

  //! the constant string type
  using const_string_t = flecsi::const_string_t;

  //! the bitfield type
  using bitfield_t = flecsi::bitfield_t;

  //! The dimension of the burton mesh.
  static constexpr size_t num_dimensions = N;

  //! The number of mesh domains in the burton mesh.
  static constexpr size_t num_domains = 2;

  //! The type for floating-point values.
  using real_t = common::real_t;

  //! The type for integer values.
  using integer_t = common::integer_t;

  //! A point type with real_t data and mesh dimension.
  using point_t = geom::point<real_t, num_dimensions>;

  //! A space ("physics") vector type with real_t data and mesh dimension.
  using vector_t = math::vector<real_t, num_dimensions>;

  //! Enumeration of the locations on the mesh where data may be sited.
  enum class attachment_site_t : size_t 
  {
    vertices,
    edges,
    faces,
    cells,
    corners,
    wedges,
    global
  };

  //! Enumeration of the available attributes on state data.
  enum class state_attribute_t : bitfield_t::field_type_t {
    persistent = 0,
    temporary = 1 // just a placeholder for now
  };

  //! \brief The locations of different bits that we set as flags
  enum bits : uint8_t
  {
    boundary = 0 //!< the boundary bit
  };

  // the tags type
  using tag_t = uint8_t;
  using tag_list_t = utils::fixed_vector< tag_t, N*N >;

  //============================================================================
  // State type definitions
  //============================================================================

  //! \brief Struct to hold meta data for private state. Meta data includes site
  //!   and attributes.
  struct private_state_meta_data_t {

    void initialize(
      attachment_site_t site_,
      bitfield_t::field_type_t attributes_) 
    {
      site = site_;
      attributes = attributes_;
    } // initialize

    attachment_site_t site;
    bitfield_t::field_type_t attributes;

  }; // struct private_state_meta_data_t

  //! A type definition of data_t based on the storage policy for the mesh.
#ifndef MESH_STORAGE_POLICY
  using data_t = flecsi::data_model::data_t<private_state_meta_data_t>;
#else
  using data_t =
    flecsi::data_model::data_t<private_state_meta_data_t, MESH_STORAGE_POLICY>;
#endif

};

} // namespace mesh
} // namespace ale

