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
 * \file tasks.h
 * 
 * \brief Simple tasks related to solving full hydro solutions.
 *
 ******************************************************************************/

#pragma once

// hydro includes
#include "types.h"

namespace apps {
namespace hydro {

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename F >
int32_t initial_conditions( T & mesh, F && ics ) {

  // get the collection accesor
  auto d = access_state( mesh, "density",    real_t );
  auto p = access_state( mesh, "pressure",   real_t );
  auto v = access_state( mesh, "velocity", vector_t );

  for ( auto c : mesh.cells() ) {
    // get the 
    auto x = c->centroid();
    // now copy the state to flexi
    std::tie( d[c], v[c], p[c] ) = std::forward<F>(ics)( x );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t update_state_from_pressure( T & mesh ) {


  // get the collection accesor
  auto eos = access_global_state( mesh, "eos", eos_t );
  state_accessor<T> state( mesh );

  for ( auto c : mesh.cells() ) {
    auto u = state(c);
    eqns_t::update_state_from_pressure( u, *eos );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t update_state_from_energy( T & mesh ) {


  // get the collection accesor
  auto eos = access_global_state( mesh, "eos", eos_t );
  state_accessor<T> state( mesh );

  for ( auto c : mesh.cells() ) {
    auto u = state(c);
    eqns_t::update_state_from_energy( u, *eos );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename T >
int32_t evaluate_time_step( T & mesh ) {

  // access what we need
  state_accessor<T> state( mesh );

  auto delta_t = access_global_state( mesh, "time_step", real_t );
  real_t cfl = access_global_state( mesh, "cfl", real_t );
 
 
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_inv(0);

  for ( auto c : mesh.cells() ) {

    // get cell properties
    auto u = state( c );
    auto vol = c->volume();

    // loop over each face
    for ( auto f : mesh.faces(c) ) {
      // estimate the length scale normal to the face
      auto delta_x = vol / f->area();
      // get the unit normal
      auto norm = f->normal();
      auto nunit = norm / abs(norm);
      // compute the inverse of the time scale
      auto dti =  E::fastest_wavespeed(u, nunit) / delta_x;
      // check for the maximum value
      dt_inv = std::max( dti, dt_inv );
    } // edge

  } // cell

  assert( dt_inv > 0 && "infinite delta t" );

  // invert dt and apply cfl
  delta_t = cfl / dt_inv;

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate residuals
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t evaluate_fluxes( T & mesh ) {

  // access what we need
  auto flux = access_state( mesh, "flux", flux_data_t );
  state_accessor<T> state( mesh );

  //----------------------------------------------------------------------------
  // TASK: loop over each edge and compute/store the flux
  // fluxes are stored on each edge
  for ( auto f : mesh.faces() ) {

    // get the normal and unit normal
    // The normal always points towards the first cell in the list.
    // So we flip it.  This makes it point from the first cell outward.
    auto norm = - f->normal();
    auto nunit = unit( norm );

    // get the cell neighbors
    auto cells = mesh.cells(f);
    auto num_cells = cells.size();


    // get the left state
    auto w_left = state( cells[0] );    

    // compute the face flux
    //
    // interior cell
    if ( num_cells == 2 ) {
      auto w_right = state( cells[1] );
      flux[f] = flux_function( w_left, w_right, nunit );
    } 
    // boundary cell
    else {
      flux[f] = boundary_flux( w_left, nunit );
    }
    
    // scale the flux by the face area
    math::multiplies_equal( flux[f], f->area() );
    
  } // for
  //----------------------------------------------------------------------------

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t apply_update( T & mesh ) {

  // access what we need
  auto flux = access_state( mesh, "flux", flux_data_t );
  state_accessor<T> state( mesh );

  // read only access
  real_t delta_t = access_global_state( mesh, "time_step", real_t );
 
  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto c : mesh.cells() ) {
    
    flux_data_t delta_u;
    math::fill( delta_u, 0.0 );

    // loop over each connected edge
    for ( auto f : mesh.faces(c) ) {
      
      // get the cell neighbors
      auto neigh = mesh.cells(f);
      auto num_neigh = neigh.size();

      // add the contribution to this cell only
      if ( neigh[0] == c )
        math::minus_equal( delta_u, flux[f] );
      else
        math::plus_equal( delta_u, flux[f] );

    } // edge

    // now compute the final update
    math::multiplies_equal( delta_u, delta_t/c->volume() );
    
    // apply the update
    auto u = state( c );
    eqns_t::update_state_from_flux( u, delta_u );

  } // for
  //----------------------------------------------------------------------------

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief Output the solution
//!
//! \param [in] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t output( T & mesh, 
                const std::string & prefix, 
                const std::string & postfix, 
                size_t output_freq ) 
{

  auto cnt = mesh.time_step_counter();
  if ( cnt % output_freq != 0 ) return 0;

  std::stringstream ss;
  ss << prefix;
  ss << std::setw( 7 ) << std::setfill( '0' ) << cnt++;
  ss << "."+postfix;
  
  mesh::write_mesh( ss.str(), mesh );
  
  return 0;
}

} // namespace hydro
} // namespace apps
