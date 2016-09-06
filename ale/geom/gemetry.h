/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some geometric functions.
////////////////////////////////////////////////////////////////////////////////

#pragma once


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

//! \brief general centroid operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for pointers
template< class T >
constexpr
decltype( std::declval<T>()->centroid() ) // exploit SFINAE
centroid( const T * t )
{ 
  return t->centroid();
}


//! \brief general centroid operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for non-pointer types
template< class T >
constexpr
decltype( std::declval<T>()->centroid() ) // exploit SFINAE
centroid( T && t )
{ 
  return std::forward<T>(t)->centroid();
}


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

//! \brief general midpoint operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for pointers
template< class T >
constexpr
decltype( std::declval<T>()->midpoint() ) // exploit SFINAE
midpoint( const T * t )
{ 
  return t->midpoint();
}


//! \brief general midpoint operator.
//! \remark this one is used for types with a centroid function
//! \remark this one is for non-pointer types
template< class T >
constexpr
decltype( std::declval<T>()->midpoint() ) // exploit SFINAE
midpoint( T && t )
{ 
  return std::forward<T>(t)->midpoint();
}

} // namespace geom
} // namespace ale
