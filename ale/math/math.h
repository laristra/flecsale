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
 * \file operators.h
 * 
 * \brief Provides a default operators for fundamental types.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <algorithm> 

// user includes
#include "../std/type_traits.h"
#include "../utils/check_types.h"
#include "detail/math.h"

namespace ale {
namespace math {

//////////////////////////////////////////////////////////////////////////////
// Compute an average
//////////////////////////////////////////////////////////////////////////////

//! \brief average operator.
//! \remark all arguments must be of the same type
//! \remark this one is used for three or more arguments
template< class T, class... Types >
constexpr
std::enable_if_t< 
  (sizeof...(Types) > 1) && 
  utils::are_type_t<T,Types...>::value, std::decay_t<T> >
average( T && t, Types&&... args )
{ 
  auto res( t ); // first one
  detail::average(res, std::forward<Types>(args)...); // sum the rest
  res /=  ( sizeof...(args) + 1 ); // divide by number
  return res;
}

//! \brief average operator.
//! \remark all arguments must be of the same type
//! \remark this one is used for two arguments only
template< class T, class U >
constexpr 
std::enable_if_t< utils::are_type_t<T,U>::value, std::decay_t<T> >
average( T && t, U && u )
{ 
  return ( std::forward<T>(t) + std::forward<T>(u) ) / 2;
}

//! \brief average operator.
//! \remark all arguments must be of the same type
//! \remark this function gives an error because not all values are the same type
template< class T, class... Types >
constexpr
std::enable_if_t< 
  !utils::are_type_t<T,Types...>::value, std::decay_t<T> >
average( T && t, Types&&... args )
{ 
  static_assert(     
    utils::are_type_t<T,Types...>::value, 
    "All parameters parameter need to be the same type." );
  return T();
}

//////////////////////////////////////////////////////////////////////////////
// return the max and min value of lists
//////////////////////////////////////////////////////////////////////////////

//! \brief return the minimum value of a list
//! \param [in] a the array to search
//! \remark general version
template< template<typename...> typename C, typename...Args >
auto min_element( const C<Args...> & a ) 
{
  return std::min_element( a.begin(), a.end() );
}

template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
auto min_element( const A<T,N...> & a ) 
{
  return std::min_element( a.begin(), a.end() );
}

//! \brief return the maximum value of a list
//! \param [in] a the array to search
//! \remark general version
template< template<typename...> typename C, typename...Args >
auto max_element( const C<Args...> & a ) 
{
  return std::max_element( a.begin(), a.end() );
}

template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
auto max_element( const A<T,N...> & a ) 
{
  return std::max_element( a.begin(), a.end() );
}

//////////////////////////////////////////////////////////////////////////////
// A dot product function
//////////////////////////////////////////////////////////////////////////////

//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template< class InputIt1, class InputIt2 >
auto dot_product( InputIt1 first1, InputIt1 last1, InputIt2 first2 )
{
  std::decay_t< decltype(*first1) > zero = 0;
  return std::inner_product(first1, last1, first2, zero );
}

//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
T dot_product(const A<T, N...> &a, const A<T, N...> &b) 
{
  auto dot = dot_product( a.begin(), a.end(), b.begin() );
  return dot;
}

template< template<typename...> typename C, typename T, typename...Args >
T dot_product(const C<T,Args...> &a, const C<T,Args...> &b) 
{
  auto dot = dot_product( a.begin(), a.end(), b.begin() );
  return dot;
}

//////////////////////////////////////////////////////////////////////////////
// magnitude of vectors
//////////////////////////////////////////////////////////////////////////////

//! \brief Compute the magnitude of the vector
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template< template<typename...> typename C, typename T, typename...Args >
T magnitude(const C<T,Args...> &a) 
{
  return std::sqrt( dot_product(a,a) );
}

template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
T magnitude(const A<T, N...> &a) 
{
  return std::sqrt( dot_product(a,a) );
}

//! \brief Compute the magnitude of the vector
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
T abs(const A<T, N...> &a) 
{
  return std::sqrt( dot_product(a,a) );
}

template< template<typename...> typename C, typename T, typename...Args >
T abs(const C<T,Args...> &a) 
{
  return std::sqrt( dot_product(a,a) );
}

//////////////////////////////////////////////////////////////////////////////
// Elementwise min and max
//////////////////////////////////////////////////////////////////////////////

//! \brief Compute the elementwise min
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
auto min(const A<T, N...> &a, const A<T, N...> &b) 
{
  A<T, N...> tmp;
  for ( auto i=0; i<tmp.size(); i++ )
    tmp[i] = std::min( a[i], b[i] );
  return tmp;
}

template< template<typename...> typename C, typename T, typename...Args >
auto min(const C<T,Args...> &a, const C<T,Args...> &b) 
{
  C<T,Args...> tmp;
  for ( auto i=0; i<a.size(); i++ )
    tmp[i] = std::min( a[i], b[i] );
  return tmp;
}


//! \brief Compute the elementwise max
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template< 
  typename T, std::size_t... N,
  template< typename, std::size_t... > typename A
 >
auto max(const A<T, N...> &a, const A<T, N...> &b) 
{
  A<T, N...> tmp;
  for ( auto i=0; i<tmp.size(); i++ )
    tmp[i] = std::max( a[i], b[i] );
  return tmp;
}

template< template<typename...> typename C, typename T, typename...Args >
auto max(const C<T,Args...> &a, const C<T,Args...> &b) 
{
  C<T,Args...> tmp;
  for ( auto i=0; i<a.size(); i++ )
    tmp[i] = std::max( a[i], b[i] );
  return tmp;
}


//////////////////////////////////////////////////////////////////////////////
// Some general math functions
//////////////////////////////////////////////////////////////////////////////

//! \brief square operator.
template< class T >
constexpr auto sqr( T && x )
{ return std::forward<T>(x) * std::forward<T>(x); }

//! \brief returns 1 if +ve, -1 if -ve
template <typename T> 
constexpr int sgn( const T & val ) {
  return (T(0) < val) - (val < T(0));
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
