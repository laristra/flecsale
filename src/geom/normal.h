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
 * \file normal.h
 * 
 * \brief Evaluate geometric normals.
 *
 ******************************************************************************/
#pragma once


// user includes
#include "ale/std/type_traits.h"
#include "ale/utils/check_types.h"


namespace ale {
namespace geom {



//! \brief general normal operator.
//! \remark this one is used for types with a normal function
//! \remark this one is for pointers
template< class T >
constexpr
decltype( std::declval<T>()->normal() ) // exploit SFINAE
normal( const T * t )
{ 
  return t->normal();
}


//! \brief general normal operator.
//! \remark this one is used for types with a normal function
//! \remark this one is for non-pointer types
template< class T >
constexpr
decltype( std::declval<T>()->normal() ) // exploit SFINAE
normal( T && t )
{ 
  return std::forward<T>(t)->normal();
}


//! \brief compute normal between two points in 2d
//! \remark all arguments must be of the same type
template< class T, class U >
constexpr
std::enable_if_t< 
  (utils::are_type_t<T,U>::value && std::decay_t<T>::size() == 2), 
  std::decay_t<T> >
normal( T && a, U && b )
{ 
  return std::decay_t<T>( a[1] - b[1], b[0] - a[0] );
}



//! \brief compute normal between two points in 3d
//! \remark all arguments must be of the same type
template< class T, class U >
constexpr
std::enable_if_t< 
  (utils::are_type_t<T,U>::value && std::decay_t<T>::size() == 3), 
  std::decay_t<T> >
normal( T && a, U && b )
{ 
  return std::decay_t<T>( a[1] * b[2] - a[2] * b[1], 
                          a[2] * b[0] - a[0] * b[2], 
                          a[0] * b[1] - a[1] * b[0] );
}

} // namespace geom
} // namespace ale
