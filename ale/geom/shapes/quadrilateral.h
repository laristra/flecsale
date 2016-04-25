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
#include "ale/math/math.h"

namespace ale {
namespace geom {

////////////////////////////////////////////////////////////////////////////////
//! the quadrilateral class
////////////////////////////////////////////////////////////////////////////////
struct quadrilateral {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::quadrilateral;

  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename... Args >
  static auto centroid( Args&&... pts ) 
  {
    return math::average( std::forward<Args>(pts)... );
  }
  
  //============================================================================
  //! \brief the volume function
  //============================================================================
  template< typename... Args >
  static auto area( Args&&... pts ) 
  {
    // auto u = pt1 - pt0;
    // auto v = pt2 - pt0;
    // auto cross = cross_product( u, v );
    // return std::abs( cross ) / 2;
    return 0.0;
  }
  
    


};

} // namespace geom
} // namespace ale
