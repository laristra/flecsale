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
 * \file tetrahedron.h
 * 
 * \brief Some utility functions for tetrahedrons.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/math/math.h"

namespace ale {
namespace geom {

////////////////////////////////////////////////////////////////////////////////
//! the tetrahedron class
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
  //! \brief the volume function
  //============================================================================
  template< typename... Args >
  static auto volume( Args&&... pts ) 
  {
    return 0.0;
  }
  
    


};

} // namespace geom
} // namespace ale
