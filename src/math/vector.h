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
//! \brief Compute the unit vector given a vector
//! \tparam T  The base value type.
//! \tparam D  The matrix/array dimension.
//! \param[in] mat  The matrix
//! \param[in] vec  The vector
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T, std::size_t D,
  template<typename, std::size_t> typename C
>
C<T,D> unit( const C<T,D> & x )
{
  auto l = abs(x);
  auto u = x / l;
  return u;
}
  

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T,
  template<typename, std::size_t> typename C
>
auto cross_product(const C<T, 2> &a, const C<T, 2> &b) 
{
  return a[0]*b[1] - a[1]*b[0];
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
