/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Some general use flux functions.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "math/general.h" 

namespace ale {
namespace eqns {


////////////////////////////////////////////////////////////////////////////////
//! \brief An averaging flux function.
//!
//! \tparam E  The equations type.
//! \tparam U  The state type.
//! \tparam V  The vector type.
//!
//! \param [in] wl,wr  The left and right states.
//! \param [in] n      The normal direction.
//! \return the flux
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
auto average_flux( const U & wl, const U & wr, const V & n) { 
  auto fl = E::flux( wl, n );
  auto fr = E::flux( wr, n );
  return fl + fr;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the rusanov flux function.
//!
//! \tparam E  the equations type
//! \tparam U  the state type
//! \tparam V  the vector type
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
//! \brief Compute the HLLE flux function.
//!
//! \tparam E  the equations type
//! \tparam U  the state type
//! \tparam V  the vector type
//!
//! \param [in] wl,wr  the left and right states
//! \param [in] n      the normal direction
//! \return the flux
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
auto hlle_flux( const U & wl, const U & wr, const V & n) { 
  // compute some things for the dissipation term
  auto sl = E::minmax_eigenvalues( wl, n );
  auto sr = E::minmax_eigenvalues( wr, n );

  auto lambda_l = std::min( sl.first, sr.first );
  auto lambda_r = std::max( sl.second, sr.second );

  if ( lambda_l >= 0 )   
    return E::flux( wl, n );
  else if ( lambda_r <= 0 ) 
    return E::flux( wr, n );
  else {
    auto fl = E::flux( wl, n );
    auto fr = E::flux( wr, n );
    auto c1 = lambda_l * lambda_r;
    auto c2 = lambda_r - lambda_l;
    auto c2inv = 1 / c2;
    auto du = E::solution_delta( wl, wr );
    //f = ( lambda_r*fl - lambda_l*fr + c1*(ur - ul) ) / c2
    decltype(du) f;
    constexpr auto num_var = f.size();
    for ( int i=0; i<num_var; ++i )
      f[i] = c2inv * ( lambda_r*fl[i] - lambda_l*fr[i] + c1*du[i] );
    return f;
  }
};


} // namespace
} // namespace

