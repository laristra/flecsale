/*~-------------------------------------------------------------------------~~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
/*!
 *
 * \file vector.h
 * 
 * \brief Provides a dimensioned array which functions as a vector.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/math/array.h"

namespace ale {
namespace math {

////////////////////////////////////////////////////////////////////////////////
//!  \brief The dimensioned_array type provides a general base for defining
//!  contiguous array types that have a specific dimension.
//!
//!  \tparam T The type of the array, e.g., P.O.D. type.
//!  \tparam D The dimension of the array, i.e., the number of elements
//!    to be stored in the array.
////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t D> 
using vector = array<T,D>;


////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the rotation matrix
//! \tparam T  The base value type.
//! \tparam D  The matrix/array dimension.
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T, std::size_t D,
  template<typename, std::size_t> typename C
>
C<T,D> reflect( const C<T,D> & v, const C<T,D> & n ) {

  auto dot = dot_product( v, n );

  auto rot = n;
  rot *= 2 * dot;

  auto tmp = v;
  tmp -= rot;

  return tmp;
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
