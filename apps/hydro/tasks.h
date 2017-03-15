/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tasks related to solving full hydro solutions.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// hydro includes
#include "types.h"

// system includes
#include <iomanip>

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
int initial_conditions( T & mesh, F && ics ) {

  // type aliases
  using counter_t = typename T::counter_t;
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;

  // get the current time
  auto soln_time = mesh.time();

  // get the collection accesor
  auto d = flecsi_get_accessor( mesh, hydro, density,    real_t, dense, 0 );
  auto p = flecsi_get_accessor( mesh, hydro, pressure,   real_t, dense, 0 );
  auto v = flecsi_get_accessor( mesh, hydro, velocity, vector_t, dense, 0 );

  auto xc = flecsi_get_accessor( mesh, mesh, cell_centroid, vector_t, dense, 0 );

  auto cs = mesh.cells();
  auto num_cells = cs.size();

  // This doesn't work with lua input
  //#pragma omp parallel for
  for ( counter_t i=0; i<num_cells; i++ ) {
    auto c = cs[i];
    std::tie( d[c], v[c], p[c] ) = std::forward<F>(ics)( xc[c], soln_time );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for updating the state using pressure.
//!
//! Updates the state from density and pressure and computes the new energy.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename EOS >
int update_state_from_pressure( T & mesh, const EOS * eos ) 
{

  // type aliases
  using counter_t = typename T::counter_t;
  using eqns_t = eqns_t<T::num_dimensions>;

  // get the collection accesor
  state_accessor<T> state( mesh );

  auto cs = mesh.cells();
  auto num_cells = cs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_cells; i++ ) {
    auto c = cs[i];
    auto u = state(c);
    eqns_t::update_state_from_pressure( u, *eos );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for updating the state from energy.
//!
//! Updates the state from density and energy and computes the new pressure.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename EOS >
int update_state_from_energy( T & mesh, const EOS * eos ) 
{

  // type aliases
  using counter_t = typename T::counter_t;
  using eqns_t = eqns_t<T::num_dimensions>;

  // get the collection accesor
  state_accessor<T> state( mesh );

  // get the cells
  auto cs = mesh.cells();
  auto num_cells = cs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_cells; i++ ) {
    auto c = cs[i];
    auto u = state(c);
    eqns_t::update_state_from_energy( u, *eos );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size.
//!
//! \tparam E  The equation of state object to use.
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename T >
int evaluate_time_step( T & mesh ) {

  // type aliases
  using  counter_t = typename T::counter_t;
  using   real_t = typename T::real_t;
  using vector_t = typename T::vector_t;

  // access what we need
  state_accessor<T> state( mesh );

  auto delta_t = flecsi_get_accessor( mesh, hydro, time_step, real_t, global, 0 );
  const auto cfl = flecsi_get_accessor( mesh, hydro, cfl, real_t, global, 0 );
 
  auto area   = mesh.face_areas();
  auto normal = mesh.face_normals();
  auto volume = mesh.cell_volumes();
 
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_inv(0);

  // get the cells
  auto cs = mesh.cells();
  auto num_cells = cs.size();

  #pragma omp parallel for reduction(max:dt_inv)
  for ( counter_t i=0; i<num_cells; i++ ) {
    auto c = cs[i];

    // get cell properties
    auto u = state( c );

    // loop over each face
    for ( auto f : mesh.faces(c) ) {
      // estimate the length scale normal to the face
      auto delta_x = volume[c] / area[f];
      // compute the inverse of the time scale
      auto dti =  E::fastest_wavespeed(u, normal[f]) / delta_x;
      // std::cout << dti << " " << delta_x << std::endl;
      // check for the maximum value
      dt_inv = std::max( dti, dt_inv );
    } // edge

  } // cell

  assert( dt_inv > 0 && "infinite delta t" );

  // invert dt and apply cfl
  *delta_t = static_cast<real_t>(cfl) / dt_inv;

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate fluxes at each face.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int evaluate_fluxes( T & mesh ) {

  // type aliases
  using counter_t = typename T::counter_t;
  using vector_t = typename T::vector_t;
  using eqns_t = eqns_t<T::num_dimensions>;
  using flux_data_t = flux_data_t<T::num_dimensions>;

  // access what we need
  auto flux = flecsi_get_accessor( mesh, hydro, flux, flux_data_t, dense, 0 );
  state_accessor<T> state( mesh );

  auto area   = mesh.face_areas();
  auto normal = mesh.face_normals();

  //----------------------------------------------------------------------------
  // TASK: loop over each edge and compute/store the flux
  // fluxes are stored on each edge
 
  // get the faces
  auto fs = mesh.faces();
  auto num_faces = fs.size();

  //for ( auto fit = fs.begin(); fit < fs.end(); ++fit  ) {

  #pragma omp parallel for
  for ( counter_t i=0; i<num_faces; i++ ) {
    auto f = fs[i];

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
      flux[f] = flux_function<eqns_t>( w_left, w_right, normal[f] );
    } 
    // boundary cell
    else {
      flux[f] = boundary_flux<eqns_t>( w_left, normal[f] );
    }
    
    // scale the flux by the face area
    flux[f] *= area[f];

    // std::cout << flux[f] << std::endl;
    
  } // for
  //----------------------------------------------------------------------------

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution in each cell.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int apply_update( T & mesh ) {

  // type aliases
  using counter_t = typename T::counter_t;
  using real_t = typename T::real_t;
  using flux_data_t = flux_data_t<T::num_dimensions>;
  using eqns_t = eqns_t<T::num_dimensions>;

  // access what we need
  auto flux = flecsi_get_accessor( mesh, hydro, flux, flux_data_t, dense, 0 );
  state_accessor<T> state( mesh );
  
  auto volume = mesh.cell_volumes();

  // read only access
  const auto delta_t = flecsi_get_accessor( mesh, hydro, time_step, real_t, global, 0 );

  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell

  auto cs = mesh.cells();
  auto num_cells = cs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_cells; i++ ) {
    
    auto c = cs[i];
    flux_data_t delta_u( 0 );

    // loop over each connected edge
    for ( auto f : mesh.faces(c) ) {
      
      // get the cell neighbors
      auto neigh = mesh.cells(f);
      auto num_neigh = neigh.size();

      // add the contribution to this cell only
      if ( neigh[0] == c )
        delta_u -= flux[f];
      else
        delta_u += flux[f];

    } // edge

    // now compute the final update
    delta_u *= static_cast<real_t>(delta_t)/volume[c];

    // apply the update
    auto u = state( c );
    eqns_t::update_state_from_flux( u, delta_u );


  } // for
  //----------------------------------------------------------------------------

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief Output the solution.
//!
//! \param [in] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int output( T & mesh, 
                const std::string & prefix, 
                const std::string & postfix, 
                size_t output_freq ) 
{

  if ( output_freq < 1 ) return 0;

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
