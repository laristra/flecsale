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
 * \file quadrilateral.h
 * 
 * \brief Some utility functions for quadrilaterals.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/geom/shapes/polygon.h"

namespace ale {
namespace geom {

////////////////////////////////////////////////////////////////////////////////
//! \brief the quadrilateral class
//! \tparam N the number of dimensions
//! \remark This is the primary template
////////////////////////////////////////////////////////////////////////////////
template<std::size_t N>
struct quadrilateral {};

////////////////////////////////////////////////////////////////////////////////
//! the quadrilateral class
//! \remark this is the 2D specialization
////////////////////////////////////////////////////////////////////////////////
template<>
struct quadrilateral<2> {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::quadrilateral;

  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename... Args >
  static auto centroid( Args&&... pts ) 
  {
    return polygon<2>::centroid( std::forward<Args>(pts)... );
  }
  
  //============================================================================
  //! \brief the volume function
  //============================================================================
  template< typename... Args >
  static auto area( Args&&... pts ) 
  { 
    return polygon<2>::area( std::forward<Args>(pts)... );
  }

  
  //============================================================================
  //! \brief the normal of the quad
  //============================================================================
  template< typename... Args >
  static auto normal( Args&&... pts ) 
  { 
    return polygon<2>::normal( std::forward<Args>(pts)... );
  }

#if 0
  
  //============================================================================
  //! \brief the normal of the quad
  //============================================================================
  template< typename T >
  static auto normal( const T & pt0, const T & pt1, const T & pt2, const T & pt3 ) 
  { 
    auto n = triangle<3>::normal( pt0, pt1, pt2 );
    n += triangle<3>::normal( pt2, pt3, pt0 );
    return n;
  }

#endif
    
};



////////////////////////////////////////////////////////////////////////////////
//! the quadrilateral class
//! \remark this is the 2D specialization
////////////////////////////////////////////////////////////////////////////////
template<>
struct quadrilateral<3> {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::quadrilateral;

  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename... Args >
  static auto centroid( Args&&... pts ) 
  {
    return polygon<3>::centroid( std::forward<Args>(pts)... );
  }


  //============================================================================
  //! \brief the midpoint function
  //============================================================================
  template< typename... Args >
  static auto midpoint( Args&&... pts ) 
  {
    return polygon<3>::midpoint( std::forward<Args>(pts)... );
  }

  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename... Args >
  static auto area( Args&&... pts ) 
  {
    return polygon<3>::area( std::forward<Args>(pts)... );
  }

  //============================================================================
  //! \brief the normal of the quad
  //============================================================================
  template< typename... Args >
  static auto normal( Args&&... pts ) 
  {
    return polygon<3>::normal( std::forward<Args>(pts)... );
  }
  

#if 0

  //============================================================================
  //! \brief the volume function
  //============================================================================
  template< typename T >
  static auto area( const T & pt0, const T & pt1, const T & pt2, const T & pt3 ) 
  { 
    auto xc = centroid( pt0, pt1, pt2, pt3 );
    auto a = 
      triangle<3>::area( pt0, pt1, xc ) + 
      triangle<3>::area( pt1, pt2, xc ) + 
      triangle<3>::area( pt2, pt3, xc ) +
      triangle<3>::area( pt3, pt0, xc );
    return a;
  }

  
  //============================================================================
  //! \brief the normal of the quad
  //============================================================================
  template< typename T >
  static auto normal( const T & pt0, const T & pt1, const T & pt2, const T & pt3 ) 
  { 
    auto xc = centroid( pt0, pt1, pt2, pt3 );
    auto n = triangle<3>::normal( pt0, pt1, xc );
    n += triangle<3>::normal( pt1, pt2, xc );
    n += triangle<3>::normal( pt2, pt3, xc );
    n += triangle<3>::normal( pt3, pt0, xc );
    return n;
  }

#endif
    
};


} // namespace geom
} // namespace ale
