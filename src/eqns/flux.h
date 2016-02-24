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
 * \file flux.h
 * 
 * \brief Some general use flux functions.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/math/math.h" 

namespace ale {
namespace eqns {


////////////////////////////////////////////////////////////////////////////////
//! \brief An averaging flux function
//!
//! \tparam E  the equations type
//!
//! \param [in] wl,wr  the left and right states
//! \param [in] n      the normal direction
//! \return the flux
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
auto average_flux( const U & wl, const U & wr, const V & n) { 
  auto fl = E::flux( wl, n );
  auto fr = E::flux( wr, n );
  return fl + fr;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the rusanov flux function
//!
//! \tparam E  the equations type
//!
//! \param [in] wl,wr  the left and right states
//! \param [in] n      the normal direction
//! \return the flux
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
auto rusanov_flux( const U & wl, const U & wr, const V & n) { 
  // get the centered flux 
  auto favg = E::flux( wl, n ) + E::flux( wr, n );
  // compute some things for the dissipation term
  auto sl = E::fastest_wavespeed( wl, n );
  auto sr = E::fastest_wavespeed( wr, n );
  auto s = std::max( sl, sr );
  auto du = E::solution_delta( wl, wr );
  // compute final flux
  // f = 0.5*(fl+fr) - s_max/2 * (ur-ul)
  math::divides_equal( favg, 2 );
  math::multiplies_equal( du, s/2 );
  return favg - du;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the HLLE flux function
//!
//! \tparam E  the equations type
//!
//! \param [in] wl,wr  the left and right states
//! \param [in] n      the normal direction
//! \return the flux
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
auto hlle_flux( const U & wl, const U & wr, const V & n) { 
  // compute some things for the dissipation term
  auto sl = E::eigenvalues( wl, n );
  auto sr = E::eigenvalues( wr, n );

  auto lambda_l = std::min( math::min_value(sl), math::min_value(sr) );
  auto lambda_r = std::max( math::max_value(sl), math::max_value(sr) );

  if ( lambda_l >= 0 )    
    return E::flux( wl, n );
  else if ( lambda_r <= 0 ) 
    return E::flux( wr, n );
  else {
    auto fl = E::flux( wl, n );
    auto fr = E::flux( wr, n );
    auto c1 = lambda_l * lambda_r;
    auto c2 = lambda_r - lambda_l;
    auto du = E::solution_delta( wl, wr );
    //f = ( lambda_r*fl - lambda_l*fr + c1*(ur - ul) ) / c2
    math::multiplies_equal( fl, lambda_r/c2 );
    math::multiplies_equal( fr, lambda_l/c2 );
    math::multiplies_equal( du, c1/c2 );
    return fl - fr + du;
  }
};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
