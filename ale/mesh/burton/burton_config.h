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
struct burton_config_t {

  //! the size type
  using size_t = std::size_t;

  //! The dimension of the burton mesh.
  static constexpr size_t num_dimensions = N;

  //! The number of mesh domains in the burton mesh.
  static constexpr size_t num_domains = 2;

  //! the constant string type
  using const_string_t = flecsi::utils::const_string_t;

  //! the bitfield type
  using bitfield_t = flecsi::utils::bitfield_t;

  //! A type used for loop indexing.
  using counter_t = long long;

  //! The type for floating-point values.
  using real_t = common::real_t;

  //! The type for integer values.
  using integer_t = common::integer_t;

  //! A point type with real_t data and mesh dimension.
  using point_t = geom::point<real_t, num_dimensions>;

  //! A space ("physics") vector type with real_t data and mesh dimension.
  using vector_t = math::vector<real_t, num_dimensions>;

  //! \brief The locations of different bits that we set as flags
  enum bits : uint8_t
  {
    boundary = 0 //!< the boundary bit
  };

  // the tags type
  using tag_t = uint8_t;
  using tag_list_t = utils::fixed_vector< tag_t, N*N >;



};

} // namespace mesh
} // namespace ale

