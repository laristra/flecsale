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
 * \file matrix.h
 * 
 * \brief Provides a dimensioned array which functions as a matrix.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/math/multi_array.h"
#include "ale/math/vector.h"

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
template <typename T, std::size_t D1, std::size_t D2> 
using matrix = multi_array<T, D1, D2>;

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T, std::size_t D,
  template<typename, std::size_t> typename C
>
auto outer_product(const C<T, D> &a, const C<T, D> &b)
{
  matrix<T,D,D> tmp;
  
  // the result is symmetric, so use the iterator to make sure we are always
  // looping in favorable order
  auto it = tmp.begin();
  
  // this order does not matter
  for ( auto i = 0; i<D; i++ )
    for ( auto j = 0; j<D; j++ )
      *it++ = a[i] * b[j];

  return tmp;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the inverse of a square matrix
//! \tparam T  The base value type.
//! \tparam D  The matrix dimension.
//! \param[in] mat  The matrix to invert
//! \return The result of the operation
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T, std::size_t D
>
auto inverse( const matrix<T, D, D> & mat )
{
  matrix<T,D,D> tmp;
  auto a = mat(0,0);
  auto b = mat(0,1);
  auto c = mat(1,0);
  auto d = mat(1,1);
  tmp(0,0) =  d;
  tmp(0,1) = -b;
  tmp(1,0) = -c;
  tmp(1,1) =  a;
  auto denom = a*d - b*c;
  assert( denom != T() );
  tmp /= denom;
  return tmp;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the product of a matrix times a vector
//! \tparam T  The base value type.
//! \tparam D  The matrix/array dimension.
//! \param[in] mat  The matrix
//! \param[in] vec  The vector
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T, std::size_t D,
  template<typename, std::size_t> typename C
>
void ax_plus_y( const matrix<T, D, D> & A, const C<T,D> & x, C<T,D> & y )
{
  for ( auto i = 0; i<D; i++ ) {
    for ( auto j = 0; j<D; j++ )
      y[i] += A(i,j) * x[j];
  }
}


////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the product of a matrix times a vector
//! \tparam T  The base value type.
//! \tparam D  The matrix/array dimension.
//! \param[in] mat  The matrix
//! \param[in] vec  The vector
////////////////////////////////////////////////////////////////////////////////
template < 
  typename T, std::size_t D,
  template<typename, std::size_t> typename C
>
auto solve( const matrix<T, D, D> & A, const C<T,D> & b )
{
  C<T,D> x(0);
  auto inv = inverse(A);
  ax_plus_y( inv, b, x );
  return x;
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
