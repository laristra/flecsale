/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utility functions for hexahedrons.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/geom/shapes/polyhedron.h"

namespace ale {
namespace geom {
namespace shapes {

////////////////////////////////////////////////////////////////////////////////
//! \brief the hexahedron class
////////////////////////////////////////////////////////////////////////////////
struct hexahedron {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::hexahedron;
  
  //============================================================================
  //! \brief the centroid function
  //============================================================================
  template< typename T >
  static auto centroid( 
    const T & pt0, const T & pt1, const T & pt2, const T & pt3,
    const T & pt4, const T & pt5, const T & pt6, const T & pt7 ) 
  {
    using face_t = std::array<T, 4>;
    polyhedron<T> poly;    
    poly.insert( face_t{pt0, pt1, pt2, pt3} );
    poly.insert( face_t{pt4, pt7, pt6, pt5} );
    poly.insert( face_t{pt0, pt4, pt5, pt1} );
    poly.insert( face_t{pt1, pt5, pt6, pt2} );
    poly.insert( face_t{pt2, pt6, pt7, pt3} );
    poly.insert( face_t{pt3, pt7, pt4, pt0} );
    return poly.centroid();
  }

  //============================================================================
  //! \brief the midpoint function
  //============================================================================
  template< typename T, typename... Args >
  static auto midpoint( T&& pt0, Args&&... pts ) 
  {
    using point_type = std::decay_t< T >;
    return polyhedron<point_type>::midpoint( 
      std::forward<T>(pt0), std::forward<Args>(pts)... );
  }
  
  //============================================================================
  //! \brief the volume function
  //! \see Grandy, "Efficient Computation of Volume of Hexahedral Cells", 1997
  //============================================================================
  template< typename T >
  static auto volume( 
    const T & pt0, const T & pt1, const T & pt2, const T & pt3,
    const T & pt4, const T & pt5, const T & pt6, const T & pt7 ) 
  {
    // broken into tets that all share the midpoint, but faster (supposedly)
    // surface is broken into tets
    auto d20 = pt2 - pt0;
    auto d50 = pt5 - pt0;
    auto d61 = pt6 - pt1;
    auto d63 = pt6 - pt3;
    auto d64 = pt6 - pt4;
    auto d70 = pt7 - pt0;
    auto det = 
      std::abs( triple_product( d61+d70, d63, d20 ) ) +
      std::abs( triple_product( d70, d63+d50, d64 ) ) +
      std::abs( triple_product( d61, d50, d64+d20 ) );
    return det / 12;
  }
     
};

} // namespace shapes
} // namespace geom
} // namespace ale
