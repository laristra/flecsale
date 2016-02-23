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
 * \file area.h
 * 
 * \brief Evaluate the area of geometric entities.
 *
 ******************************************************************************/
#pragma once

// user includes
#include "ale/std/type_traits.h"
#include "ale/utils/check_types.h"
#include "detail/area.h"

namespace ale {
namespace geom {

//! \brief general area operator.
//! \remark this one is used for types with an area function
//! \remark this one is for pointers
template< class T >
constexpr
decltype( std::declval<T>()->area() ) // exploit SFINAE
area( const T * t )
{ 
  return t->area();
}


//! \brief general area operator.
//! \remark this one is used for types with a area function
//! \remark this one is for non-pointer types
template< class T >
constexpr
decltype( std::declval<T>()->area() ) // exploit SFINAE
area( T && t )
{ 
  return std::forward<T>(t)->area();
}


//! \brief compute area for 2d
//! \remark all arguments must be of the same type
//! \remark this one is used for three or more arguments
template< class T, class... Types >
constexpr
std::enable_if_t< 
  (sizeof...(Types) > 2 && utils::are_type_t<T,Types...>::value && 
   std::decay_t<T>::size() == 2), 
  std::decay_t<T> >
area( T && t, Types&&... args )
{ 
  // initialize area
  typename std::decay_t<T>::value_type area(0);
  // call main implementation, sticking first point back on end
  detail::area_2d( area, 
                   std::forward<T>(t), 
                   std::forward<Types>(args)..., 
                   std::forward<T>(t) );
  return std::abs(area/2);
}


//! \brief compute area for 2d or 3d
//! \remark this is a runtime computation
template< 
  typename T, typename... Args, 
  template <typename, typename...> class V 
>
auto area_2d( const V<T, Args...> & points )
{ 
  assert( points.size() > 2 && "not enough points for 2d" );

  T cx(0);
  typename T::value_type a(0);

  auto n = points.size();

  auto po = &points.back();
  for ( auto i=0; i<n; i++ ) {
    auto pn = &points[i];
    a += (*po)[0]*(*pn)[1] - (*pn)[0]*(*po)[1];
    po = pn;
  }
  return std::abs(a/2);
}


//! \brief compute area for 2d or 3d
//! \remark this is a runtime computation
template< 
  typename T, typename... Args, 
  template <typename, typename...> class V 
>
auto area( const V<T, Args...> & points )
{ 
  assert( points.size() > 0 && "not enough points" );
  auto dim = points[0].size();
  if ( dim == 2 )
    return area_2d( points );
  else
    assert( false && "unsuported dimemsion" );
}


} // namespace geom
} // namespace ale
