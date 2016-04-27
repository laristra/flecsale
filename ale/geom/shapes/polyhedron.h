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
 * \file polyhedron.h
 * 
 * \brief Some utility functions for polyhedrons.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/math/math.h"
#include "ale/utils/array_ref.h"

namespace ale {
namespace geom {

////////////////////////////////////////////////////////////////////////////////
//! \brief the polyhedron class
//! \see Stroud, Approximate calculation of multiple integrals, 
//!      Prentice-Hall Inc., 1971.
////////////////////////////////////////////////////////////////////////////////
template< typename P >
class polyhedron {

public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::polyhedron;

  //! \brief the point and coordinate types
  using point_type = P;
  using coord_type = typename point_type::value_type;

  //============================================================================
  //! \brief insert a face into the polyhedron
  //============================================================================
  template< typename T >
  void insert( utils::array_ref<T> points ) 
  {
    insert( points.begin(), points.end() );
  }

  template< typename T >
  void insert( std::initializer_list<T> points ) 
  {
    insert( points.begin(), points.end() );
  }

  template< typename InputIt >
  void insert( InputIt first, InputIt last ) 
  {
    // copy the points
    std::vector< point_type > points( first, last );
    // move into the vector
    faces_.emplace_back( std::move(points) );
  }


  //============================================================================
  //! \brief the volume function
  //============================================================================
  auto centroid() 
  {
    // initialize volume
    point_type cx(0);
    coord_type v = 0;

    //--------------------------------------------------------------------------
    // loop over faces
    for ( const auto & points : faces_ ) {

      // face midpoint
      auto xm = math::average( points );

      // for each face edge
      auto po = std::prev( points.end() );
      for ( auto pn=points.begin(); pn!=points.end(); pn++ ) {
        // get normal
        auto n = triangle<3>::normal( *po, *pn, xm );
        // compute main contribution
        auto a1 = *po + *pn;
        auto a2 = *pn +  xm;
        auto a3 =  xm + *po;
        a1 *= a1;
        a2 *= a2;
        a3 *= a3;
        auto prod = a1;
        prod += a2;
        prod += a3;
        // multiply by the normal
        prod *= n;
        // add contribution to centroid
        cx += prod;
        // dot with any coordinate for volume
        v += dot_product( n, xm );
        // store old point
        po = pn;
      }
        
    }

    //--------------------------------------------------------------------------
    // return result

    // divide by volume
    cx /= 8 * v;

    return cx;
  }

  //============================================================================
  //! \brief the volume function
  //============================================================================
  auto volume() 
  {

    // initialize volume
    coord_type v = 0;

    //--------------------------------------------------------------------------
    // loop over faces
    for ( const auto & points : faces_ ) {

      // face midpoint
      auto xm = math::average( points );

      // for each face edge
      auto po = std::prev( points.end() );
      for ( auto pn=points.begin(); pn!=points.end(); pn++ ) {
        // get normal
        auto n = triangle<3>::normal( *po, *pn, xm );
        // dot with any coordinate
        v += dot_product( n, xm );
        // store old point
        po = pn;
      }
        
    }

    //--------------------------------------------------------------------------
    // return result

    return std::abs(v) / 3;
  }
  
    
  //============================================================================
  // Private Data
  //============================================================================
private:

  //! the coordinates of each face
  std::vector< std::vector<point_type> > faces_;

};

} // namespace geom
} // namespace ale
