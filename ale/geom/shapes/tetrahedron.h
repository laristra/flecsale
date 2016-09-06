/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utility functions for tetrahedrons.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/math/math.h"

namespace ale {
namespace geom {
namespace shapes {

////////////////////////////////////////////////////////////////////////////////
//! \brief the tetrahedron class
////////////////////////////////////////////////////////////////////////////////
struct tetrahedron {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::tetrahedron;

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
  static auto volume( 
    const T & pt0, const T & pt1, const T & pt2, const T & pt3 ) 
  {
    auto det = 
      pt0[0]*pt1[1]*pt2[2] + 
      pt0[1]*pt1[2]*pt3[0] + 
      pt0[2]*pt2[0]*pt3[1] + 
      pt1[0]*pt2[1]*pt3[2] -
      pt3[0]*pt2[1]*pt1[2] - 
      pt3[1]*pt2[2]*pt0[0] - 
      pt3[2]*pt1[0]*pt0[1] - 
      pt2[0]*pt1[1]*pt0[2];
    return std::abs(det) / 6;
  }
  
    


};

} // namespace shapes
} // namespace geom
} // namespace ale
