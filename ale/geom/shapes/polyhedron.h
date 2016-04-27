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

namespace ale {
namespace geom {

////////////////////////////////////////////////////////////////////////////////
//! \brief the polyhedron class
//! \see Stroud, Approximate calculation of multiple integrals, 
//!      Prentice-Hall Inc., 1971.
////////////////////////////////////////////////////////////////////////////////
struct polyhedron {

  //! \brief the shape identifier
  static constexpr auto shape = geometric_shapes_t::polyhedron;

  //============================================================================
  //! \brief the volume function
  //============================================================================
  template< typename P >
  static
  auto centroid( std::initializer_list< std::initializer_list<P> > faces ) 
  {
    using point_type = typename std::decay_t<P>;
    using coord_type = typename point_type::value_type;

    // initialize volume
    point_type cx(0);
    coord_type v = 0;

    //--------------------------------------------------------------------------
    // loop over faces
    for ( const auto & points : faces ) {

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
        auto prod = a1 * a1;
        prod += a2 * a2;
        prod += a3 * a3;
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
  template<
    typename P, typename... Args, 
    template<typename,typename...> typename V
  >
  static
  auto volume( const V<P,Args...> & faces ) 
  {
    using point_type = typename std::decay_t<P>::value_type;
    using coord_type = typename point_type::value_type;

    // initialize volume
    coord_type v = 0;

    //--------------------------------------------------------------------------
    // loop over faces
    for ( const auto & points : faces ) {

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
    return std::abs(v) / 6;
  }
  
    


};

} // namespace geom
} // namespace ale
