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
#include "ale/utils/errors.h"

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
//! \brief Compute the normal between two points
//! \tparam T  The array base value type.
//! \param[in] a  The first point
//! \param[in] b  The other point
//! \return The result of the operation
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T,
  template<typename, std::size_t> typename C
>
C<T, 2> normal(const C<T, 2> &a, const C<T, 2> &b) 
{
  return { a[1] - b[1], b[0] - a[0] };
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the normal between two points
//! \tparam T  The array base value type.
//! \param[in] a  The first point
//! \param[in] b  The other point
//! \return The result of the operation
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T,
  template<typename, std::size_t> typename C
>
C<T, 3> normal(const C<T, 3> &a, const C<T, 3> &b) 
{
  raise_runtime_error("you should never get here");
  return { 0, 0, 0 }; // FIXME - this is here as a hack
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the cross product
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
T cross_product(const C<T, 2> &a, const C<T, 2> &b) 
{
  return a[0]*b[1] - a[1]*b[0];
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the cross product
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
auto cross_product(const C<T, 3> &a, const C<T, 3> &b) 
{
  C<T, 3> tmp;
  tmp[0] = a[1]*b[2] - a[2]*b[1];
  tmp[1] = a[2]*b[0] - a[0]*b[2];
  tmp[2] = a[0]*b[1] - a[1]*b[0];
  return tmp;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the cross product
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
T triple_product(const C<T, 3> &a, const C<T, 3> &b, const C<T, 3> &c) 
{
  return 
    a[0]*b[1]*c[2] + b[0]*c[1]*a[2] + c[0]*a[1]*b[2] -
    a[2]*b[1]*c[0] - b[2]*c[1]*a[0] - c[2]*a[1]*b[0];;
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
