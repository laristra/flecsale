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
  auto M = access_state( mesh, "cell:mass",       real_t );
  auto V = access_state( mesh, "cell:volume",     real_t );
  auto p = access_state( mesh, "cell:pressure",   real_t );
  auto v = access_state( mesh, "cell:velocity", vector_t );

  for ( auto c : mesh.cells() ) {
    // get the 
    auto x = c->centroid();
    // now copy the state to flexi
    real_t d;
    std::tie( d, v[c], p[c] ) = std::forward<F>(ics)( x );
    // set mass and volume now
    auto vol = c->area();
    M[c] = d*vol;
    V[c] = vol;
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
  cell_state_accessor<T> cell_state( mesh );

  for ( auto c : mesh.cells() ) {
    auto u = cell_state(c);
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
  cell_state_accessor<T> cell_state( mesh );

  for ( auto c : mesh.cells() ) {
    auto u = cell_state(c);
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
  cell_state_accessor<T> cell_state( mesh );

  auto delta_t = access_global_state( mesh, "time_step", real_t );
  real_t cfl = access_global_state( mesh, "cfl", real_t );
 
 
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_inv(0);

  for ( auto c : mesh.cells() ) {

    // get cell properties
    auto u = cell_state( c );
    auto area = c->area();

    // loop over each edge
    for ( auto e : mesh.edges(c) ) {
      // estimate the length scale normal to the edge
      auto delta_x = area / e->length();
      // get the unit normal
      auto norm = e->normal();
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
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t evaluate_nodal_state( T & mesh ) {

  // access what we need
  vertex_state_accessor<T> vertex_state( mesh );
  cell_state_accessor<T>     cell_state( mesh );

  //----------------------------------------------------------------------------
  // Loop over each vertex
  for ( auto v : mesh.vertices() ) {
    cout << v.id() << " : ";
    // loop over each connected edge
    for ( auto c : mesh.cells(v) ) {
      cout << c.id() << " ";
    } // cell
    cout << endl;
  } // vertex
  //----------------------------------------------------------------------------

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
  cell_state_accessor<T> cell_state( mesh );

  //----------------------------------------------------------------------------
  // TASK: loop over each edge and compute/store the flux
  // fluxes are stored on each edge
  for ( auto e : mesh.edges() ) {

    // get the normal and unit normal
    // The normal always points towards the first cell in the list.
    // So we flip it.  This makes it point from the first cell outward.
    auto norm = - e->normal();
    auto nunit = norm / abs(norm);

    // get the cell neighbors
    auto c = mesh.cells(e).to_vec();
    auto cells = c.size();


    // get the left state
    auto w_left = cell_state( c[0] );    

    // interior cell
    if ( cells == 2 ) {
      // Normal always points from the first cell to the second
      // dot product below should always be the same sign ( positive )
      auto cx_left = c[0]->centroid();
      auto cx_right = c[1]->centroid();
      auto delta = cx_right - cx_left;
      auto dot = math::dot_product( norm, delta );
      if ( dot < 0 ) norm = - norm;
      assert( dot >= 0 && "interior normal points wrong way!!!!" );
      // compute the interface flux
      auto w_right = cell_state( c[1] );
      flux[e] = flux_function( w_left, w_right, nunit );
    } 
    // boundary cell
    else {
      // Normal always points outside, so the result of this
      // dot product below should always be the same sign ( positive )
      auto mp = e->midpoint();
      auto cx = c[0]->centroid();
      auto delta = mp - cx;
      auto dot = math::dot_product( norm, delta );
      if ( dot < 0 ) norm = - norm;
      assert( dot >= 0 && "boundary normal points inward!!!!" );
      // compute the boundary flux
      flux[e] = boundary_flux( w_left, nunit );
    }
    
    // scale the flux by the face area
    math::multiplies_equal( flux[e], e->length() );
    
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
  cell_state_accessor<T> cell_state( mesh );

  // read only access
  real_t delta_t = access_global_state( mesh, "time_step", real_t );
 
  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto cell : mesh.cells() ) {
    
    flux_data_t delta_u;
    math::fill( delta_u, 0.0 );

    // loop over each connected edge
    for ( auto edge : mesh.edges(cell) ) {
      
      // get the cell neighbors
      auto neigh = mesh.cells(edge).to_vec();
      auto num_neigh = neigh.size();

      // add the contribution to this cell only
      if ( neigh[0] == cell )
        math::minus_equal( delta_u, flux[edge] );
      else
        math::plus_equal( delta_u, flux[edge] );

    } // edge

    // now compute the final update
    math::multiplies_equal( delta_u, delta_t/cell->area() );
    
    // apply the update
    auto u = cell_state( cell );
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
int32_t output( T & mesh, const std::string & prefix ) 
{
  static size_t cnt = 0;

  std::stringstream ss;
  ss << prefix;
  ss << std::setw( 7 ) << std::setfill( '0' ) << cnt++;
  ss << ".exo";
  
  mesh::write_mesh( ss.str(), mesh );
  
  return 0;
}

} // namespace hydro
} // namespace apps
