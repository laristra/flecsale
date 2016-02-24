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


// user includes
#include "ale/std/type_traits.h"
#include "ale/utils/check_types.h"
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
