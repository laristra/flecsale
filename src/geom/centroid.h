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
 * \file centroid.h
 * 
 * \brief Evaluate the centroid of a list of points.
 *
 ******************************************************************************/
#pragma once


// user includes
#include "ale/std/type_traits.h"
#include "ale/utils/check_types.h"
#include "detail/centroid.h"


namespace ale {
namespace geom {



//! \brief general centroid operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for pointers
template< class T >
constexpr
auto centroid( const T * t )
{ 
  return t->centroid();
}


//! \brief general centroid operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for non-pointer types
template< class T >
constexpr
auto centroid( T && t )
{ 
  return std::forward<T>(t)->centroid();
}


//! \brief compute centroid for 2d
//! \remark all arguments must be of the same type
//! \remark this one is used for three or more arguments
template< class T, class... Types >
constexpr
std::enable_if_t< 
  (sizeof...(Types) > 2) && utils::are_type_t<T,Types...>::value, 
  std::decay_t<T> >
centroid_2d( T && t, Types&&... args )
{ 
  // initialize centroid and volume
  std::decay_t<T> cx( 0 );
  typename std::decay_t<T>::value_type vol(0);
  // call main implementation, sticking first point back on end
  detail::centroid_2d( cx, vol, 
                       std::forward<T>(t), 
                       std::forward<Types>(args)..., 
                       std::forward<T>(t) );
  cx /=  6 * vol; // divide by volume
  return cx;
}




} // namespace geom
} // namespace ale
