/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file triangle.h
/// \brief Some utility functions for triangles.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/math/general.h"

namespace ale {
namespace geom {
namespace shapes {

////////////////////////////////////////////////////////////////////////////////
//! \brief the triangle class
//! \tparam N the number of dimensions
//! \remark This is the primary template
////////////////////////////////////////////////////////////////////////////////
template<std::size_t N>
struct triangle {};

////////////////////////////////////////////////////////////////////////////////
//! \brief the triangle class
//! \remark this is the 2D specialization
////////////////////////////////////////////////////////////////////////////////
template<>
struct triangle<2> {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::triangle;

  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename... Args >
  static auto centroid( Args&&... pts ) 
  {
    return math::average( std::forward<Args>(pts)... );
  }

  //============================================================================
  //! \brief the midpoint is the same as the centroid
  //============================================================================
  template< typename... Args >
  static auto midpoint( Args&&... pts ) 
  {
    return math::average( std::forward<Args>(pts)... );
  }

  //============================================================================
  //! \brief the volume function
  //============================================================================
  template< typename T >
  static auto area( const T & pt0, const T & pt1, const T & pt2 ) 
  {
    auto u = pt1 - pt0;
    auto v = pt2 - pt0;
    auto cross = cross_product( u, v );
    return std::abs( cross ) / 2;
  }
  
  //============================================================================
  //! \brief the normal of the triangle
  //============================================================================
  template< typename T >
  static auto normal( const T & pt0, const T & pt1, const T & pt2 ) 
  {
    auto u = pt1 - pt0;
    auto v = pt2 - pt0;
    auto cross = cross_product( u, v );
    cross /= 2;
    return T{0, cross};
  }
  
};



////////////////////////////////////////////////////////////////////////////////
//! \brief the triangle class
//! \remark this is the 3D specialization
////////////////////////////////////////////////////////////////////////////////
template<>
struct triangle<3> {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::triangle;

  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename... Args >
  static auto centroid( Args&&... pts ) 
  {
    return math::average( std::forward<Args>(pts)... );
  }


  //============================================================================
  //! \brief the midpoint is the same as the centroid
  //============================================================================
  template< typename... Args >
  static auto midpoint( Args&&... pts ) 
  {
    return math::average( std::forward<Args>(pts)... );
  }
  
  //============================================================================
  //! \brief the volume function
  //============================================================================
  template< typename T >
  static auto area( const T & pt0, const T & pt1, const T & pt2 ) 
  {
    auto u = pt1 - pt0;
    auto v = pt2 - pt0;
    auto cross = cross_product( u, v );
    return math::abs( cross ) / 2;
  }

  //============================================================================
  //! \brief the normal of the trianlge
  //============================================================================
  template< typename T >
  static auto normal( const T & pt0, const T & pt1, const T & pt2 ) 
  {
    auto u = pt1 - pt0;
    auto v = pt2 - pt0;
    auto cross = cross_product( u, v );
    cross /= 2;
    return cross;
  }
  
};

} // namespace shapes
} // namespace geom
} // namespace ale
