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

// system includes
#include <algorithm>
#include <array>
#include <functional> 

//! user includes
#include "ale/utils/check_types.h"

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
template <typename T, size_t D> 
using vector_t = std::array<T,D>;


//! \brief Add to a vector_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename T, size_t D>
void add_to( vector_t<T,D>& lhs, 
             const vector_t<T,D>& rhs )
{
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
                  [](const auto & a, const auto & b) { return a + b; } );
}

template <typename T, size_t D>
void add_to( vector_t<T,D>& lhs, 
             const auto& rhs )
{
  std::transform( lhs.begin(), lhs.end(), lhs.begin(),
                  [&](const auto & a) { return a + rhs; } );
}



//! \brief Addition operator involving vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename T, size_t D>
auto operator+( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), tmp.begin(), 
                  [](const auto & a, const auto & b) { return a + b; } );
  return tmp;
}

template <typename T, size_t D>
auto operator+( const vector_t<T,D>& lhs, 
                const auto& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(), 
                  [&](const auto & a) { return a + rhs; } );
  return tmp;
}

template <typename T, size_t D>
auto operator+( const auto& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(), 
                  [&](const auto & a) { return lhs + a; } );
  return tmp;
}




//! \brief Subtract from a vector_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename T, size_t D>
void subtract_from( vector_t<T,D>& lhs, 
                    const vector_t<T,D>& rhs )
{
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
                  [](const auto & a, const auto & b) { return a - b; } );
}

template <typename T, size_t D>
void subtract_from( vector_t<T,D>& lhs, 
                    const auto& rhs )
{
  std::transform( lhs.begin(), lhs.end(), lhs.begin(),
                  [&](const auto & a) { return a - rhs; } );
}



//! \brief Minus operator involving vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename T, size_t D>
auto operator-( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), tmp.begin(), 
                  [](const auto & a, const auto & b) { return a - b; } );
  return tmp;
}

template <typename T, size_t D>
auto operator-( const vector_t<T,D>& lhs, 
                const auto& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(), 
                  [&](const auto & a) { return a - rhs; } );
  return tmp;
}

template <typename T, size_t D>
auto operator-( const auto& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(), 
                  [&](const auto & a) { return lhs - a; } );
  return tmp;
}


//! \brief multiply a vector_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename T, size_t D>
void multiply_by( vector_t<T,D>& lhs, 
                  const vector_t<T,D>& rhs )
{
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
                  [](const auto & a, const auto & b) { return a * b; } );
}

template <typename T, size_t D>
void multiply_by( vector_t<T,D>& lhs, 
                  const auto& rhs )
{
  std::transform( lhs.begin(), lhs.end(), lhs.begin(),
                  [&](const auto & a) { return a * rhs; } );
}



//! \brief Multiply operator involving vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename T, size_t D>
auto operator*( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), tmp.begin(), 
                  [](const auto & a, const auto & b) { return a * b; } );
  return tmp;
}

template <typename T, size_t D>
auto operator*( const vector_t<T,D>& lhs, 
                const auto& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(), 
                  [&](const auto & a) { return a * rhs; } );
  return tmp;
}

template <typename T, size_t D>
auto operator*( const auto& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(), 
                  [&](const auto & a) { return lhs * a; } );
  return tmp;
}


//! \brief divide a vector_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename T, size_t D>
void divide_by( vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
                  [](const auto & a, const auto & b) { return a / b; } );
}

template <typename T, size_t D>
void divide_by( vector_t<T,D>& lhs, 
                const auto& rhs )
{
  std::transform( lhs.begin(), lhs.end(), lhs.begin(),
                  [&](const auto & a) { return a / rhs; } );
}



//! \brief Division operator involving vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename T, size_t D>
auto operator/( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), tmp.begin(), 
                  [](const auto & a, const auto & b) { return a / b; } );
  return tmp;
}

template <typename T, size_t D>
auto operator/( const vector_t<T,D>& lhs, 
                const auto& rhs )
{
  vector_t<T,D> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(), 
                  [&](const auto & a) { return a / rhs; } );
  return tmp;
}

template <typename T, size_t D>
auto operator/( const auto& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(), 
                  [&](const auto & a) { return lhs / a; } );
  return tmp;
}




//! \brief Output operator for vector_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in,out] os  The ostream to dump output to.
//! \param[in]     rhs The vector_t on the right hand side of the operator.
//! \return A reference to the current ostream.
template <typename T, size_t D>
auto & operator<<(std::ostream& os, const vector_t<T,D>& a)
{
  os << "(";
  for ( size_t i=0; i<D; i++ ) 
    os << " " << a[i];
  os << " )";
  return os;
}





//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename T, size_t D>
auto dot_product(const vector_t<T, D> &a, const vector_t<T, D> &b) {

  const auto zero = T();
  return std::inner_product(a.begin(), a.end(), b.begin(), zero);
}

template <typename T>
auto dot_product(const vector_t<T, 1> &a, const vector_t<T, 1> &b) {
  return a[0] * b[0];
}

//! \brief Compute the magnitude of the vector
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename T, size_t D> 
auto magnitude(const vector_t<T, D> &a) {
  return std::sqrt( dot_product(a,a) );
}

   

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
