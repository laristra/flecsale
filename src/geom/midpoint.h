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
 * \file midpoint.h
 * 
 * \brief Evaluate the midpoint of two points.
 *
 ******************************************************************************/
#pragma once


namespace ale {
namespace geom {

//! \brief general midpoint operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for pointers
template< class T >
constexpr
auto midpoint( const T * t )
{ 
  return t->midpoint();
}


//! \brief general midpoint operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for non-pointer types
template< class T >
constexpr
auto midpoint( T && t )
{ 
  return std::forward<T>(t).midpoint();
}



} // namespace geom
} // namespace ale
