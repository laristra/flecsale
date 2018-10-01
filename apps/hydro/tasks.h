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

// flecsi includes
#include <flecsale/io/io_exodus.h>
#include <flecsi/execution/context.h>
#include <flecsi/execution/execution.h>
#include <ristra/utils/string_utils.h>

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
void initial_conditions( 
  client_handle_r__<mesh_t>  mesh,
  inputs_t::ics_function_t ics, 
  eos_t eos,
  real_t soln_time,
  dense_handle_w__<real_t> d,
  dense_handle_w__<vector_t> v,
  dense_handle_w__<real_t> e,
  dense_handle_w__<real_t> p,
  dense_handle_w__<real_t> T,
  dense_handle_w__<real_t> a
) {

  // This doesn't work with lua input
  //#pragma omp parallel for
  for ( auto c : mesh.cells( flecsi::owned ) ) {
    std::tie( d(c), v(c), p(c) ) = ics( c->centroid(), soln_time );
    eqns_t::update_state_from_pressure( 
      pack( c, d, v, p, e, T, a ),
      eos
    );
  }

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size.
//!
//! \tparam E  The equation of state object to use.
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
real_t evaluate_time_step(
  client_handle_r__<mesh_t> mesh,
  dense_handle_r__<real_t> d,
  dense_handle_r__<vector_t> v,
  dense_handle_r__<real_t> e,
  dense_handle_r__<real_t> p,
  dense_handle_r__<real_t> T,
  dense_handle_r__<real_t> a,
  real_t CFL,
  real_t max_dt
) {
 
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_inv(0);

  for ( auto c : mesh.cells( flecsi::owned ) ) {

    // get the solution state
    auto u = pack( c, d, v, p, e, T, a );

    // loop over each face
    for ( auto f : mesh.faces(c) ) {
      // estimate the length scale normal to the face
      auto delta_x = c->volume() / f->area();
      // compute the inverse of the time scale
      auto dti = eqns_t::fastest_wavespeed( u, f->normal() ) / delta_x;
      // check for the maximum value
      dt_inv = std::max( dti, dt_inv );
    } // edge

  } // cell

  if ( dt_inv <= 0 ) 
    throw_runtime_error( "infinite delta t" );

  real_t time_step = 1 / dt_inv;
  time_step *= CFL;

  // access the computed time step and make sure its not too large
  time_step = std::min( time_step, max_dt );

  return time_step;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate fluxes at each face.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void evaluate_fluxes( 
  client_handle_r__<mesh_t> mesh,
  dense_handle_r__<real_t> d,
  dense_handle_r__<vector_t> v,
  dense_handle_r__<real_t> e,
  dense_handle_r__<real_t> p,
  dense_handle_r__<real_t> T,
  dense_handle_r__<real_t> a,
  dense_handle_w__<flux_data_t> flux
) {

  const auto & face_list = mesh.faces( flecsi::owned );
  auto num_faces = face_list.size();

  #pragma omp parallel for
  for ( counter_t fit = 0; fit < num_faces; ++fit )
  {

    const auto & f = face_list[fit];
    
    // get the cell neighbors
    const auto & cells = mesh.cells(f);
    auto num_cells = cells.size();

    // get the left state
    auto w_left = pack( cells[0], d, v, p, e, T, a );
    
    // compute the face flux
    //
    // interior cell
    if ( num_cells == 2 ) {
      auto w_right = pack( cells[1], d, v, p, e, T, a );
      flux(f) = flux_function<eqns_t>( w_left, w_right, f->normal() );
    } 
    // boundary cell
    else {
      flux(f) = boundary_flux<eqns_t>( w_left, f->normal() );
    }
   
    // scale the flux by the face area
    flux(f) *= f->area();

  } // for
  //----------------------------------------------------------------------------

}

////////////////////////////////////////////////////////////////////////////////
//! \brief Perform a reduction to get the time step
////////////////////////////////////////////////////////////////////////////////
void gather_time_step(
  future_handle__<real_t> local_time_step,
  global_handle_w__<real_t> time_step
) {

  auto dt = flecsi::execution::context_t::instance().reduce_min(local_time_step);
  time_step = 0.;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution in each cell.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void apply_update( 
  client_handle_r__<mesh_t> mesh,
  eos_t eos,
  real_t delta_t,
  dense_handle_r__<flux_data_t> flux,
  dense_handle_rw__<real_t> d,
  dense_handle_rw__<vector_t> v,
  dense_handle_rw__<real_t> e,
  dense_handle_rw__<real_t> p,
  dense_handle_rw__<real_t> T,
  dense_handle_rw__<real_t> a
) {

  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell


  //auto delta_t = static_cast<real_t>( time_step );

  const auto & cell_list = mesh.cells( flecsi::owned );
  auto num_cells = cell_list.size();

  #pragma omp parallel for
  for ( counter_t cit = 0; cit < num_cells; ++cit )
  {

    const auto & c = cell_list[cit];

    // initialize the update
    flux_data_t delta_u( 0 );

    // loop over each connected edge
    for ( auto f : mesh.faces(c) ) {
      
      // get the cell neighbors
      auto neigh = mesh.cells(f);
      auto num_neigh = neigh.size();

      // add the contribution to this cell only
      if ( neigh[0] == c )
        delta_u -= flux(f);
      else
        delta_u += flux(f);

    } // edge

    // now compute the final update
    delta_u *= delta_t/c->volume();

    // apply the update
    auto u = pack(c, d, v, p, e, T, a);
    eqns_t::update_state_from_flux( u, delta_u );

    // update the rest of the quantities
    eqns_t::update_state_from_energy( u, eos );

    // check the solution quantities
    if ( eqns_t::internal_energy(u) < 0 || eqns_t::density(u) < 0 ) 
      throw_runtime_error( "Negative density or internal energy encountered!" );

  } // for
  //----------------------------------------------------------------------------
}


////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
////////////////////////////////////////////////////////////////////////////////
void output( 
  client_handle_r__<mesh_t> mesh, 
  char_array_t prefix,
	char_array_t postfix,
	size_t iteration,
	real_t time,
  dense_handle_r__<real_t> d,
  dense_handle_r__<vector_t> v,
  dense_handle_r__<real_t> e,
  dense_handle_r__<real_t> p,
  dense_handle_r__<real_t> T,
  dense_handle_r__<real_t> a
) {
  clog(info) << "OUTPUT MESH TASK" << std::endl;
 
  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // figure out this ranks file name
  auto output_filename = 
    prefix.str() + "_rank" + apps::common::zero_padded(rank) +
    "." + apps::common::zero_padded(iteration) + "." + postfix.str();

  // now outut the mesh
  flecsale::io::io_exodus__<mesh_t>::write(
    output_filename, mesh, iteration, time, &d //, v, e, p, T, a
  );
}

////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
////////////////////////////////////////////////////////////////////////////////
void print( 
  client_handle_r__<mesh_t> mesh,
  char_array_t filename
) {

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  clog(info) << "PRINT MESH ON RANK " << rank << std::endl;
 
  // figure out this ranks file name
  auto name_and_ext = ristra::utils::split_extension( filename.str() );
  auto output_filename = 
    name_and_ext.first + "_rank" + apps::common::zero_padded(rank) +
    "." + name_and_ext.second;

  // dump to file
  std::cout << "Dumping connectivity to: " << output_filename << std::endl;
  std::ofstream file( output_filename );
  mesh.dump( file );

  // close file
  file.close();
  
}


////////////////////////////////////////////////////////////////////////////////
// TASK REGISTRATION
////////////////////////////////////////////////////////////////////////////////

flecsi_register_task(initial_conditions, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(evaluate_time_step, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(gather_time_step, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(evaluate_fluxes, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(apply_update, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(output, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(print, apps::hydro, loc, single|flecsi::leaf);

} // namespace hydro
} // namespace apps
