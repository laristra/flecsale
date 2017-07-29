/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This is the main driver for the hydro solver.
///////////////////////////////////////////////////////////////////////////////
#pragma once

// hydro includes
#include "arguments.h"
#include "types.h"

// user includes
//#include <flecsale/mesh/mesh_utils.h>
#include <flecsale/utils/time_utils.h>
#include <flecsale/io/catalyst/adaptor.h>


// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

namespace apps {
namespace hydro {
  
// type aliases
using flux_data_t = flux_data__<mesh_t::num_dimensions>;

// register the data client
flecsi_register_data_client(mesh_t, meshes, mesh0);

// create some field data.  Fields are registered as struct of arrays.
// this allows us to access the data in different patterns.
flecsi_register_field(
  mesh_t, 
  hydro,  
  density,   
  mesh_t::real_t, 
  dense, 
  2, 
  mesh_t::index_spaces_t::cells
);

// flecsi_register_field(mesh_t, hydro, pressure,   mesh_t::real_t, dense, 1, cells);
// flecsi_register_field(mesh_t, hydro, velocity, mesh_t::vector_t, dense, 2, cells);

// flecsi_register_field(mesh_t, hydro, internal_energy, mesh_t::real_t, dense, 2, cells);
// flecsi_register_field(mesh_t, hydro,     temperature, mesh_t::real_t, dense, 1, cells);
// flecsi_register_field(mesh_t, hydro,     sound_speed, mesh_t::real_t, dense, 1, cells);

// compute the fluxes.  here I am regestering a struct as the stored data
// type since I will only ever be accesissing all the data at once.
// flecsi_register_field(mesh_t, hydro, flux, flux_data_t, dense, 1, faces);

// register the time step and set a cfl
//flecsi_register_field( mesh_t, hydro, time_step, real_t, global, 1 );
//flecsi_register_field( mesh_t, hydro, cfl, real_t, global, 1 );

// Register the total energy
//flecsi_register_field( mesh_t, hydro, sum_total_energy, real_t, global, 1 );


///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int driver(int argc, char** argv) 
{

  //===========================================================================
  // Parse arguments
  //===========================================================================
  
  auto args = process_arguments( argc, argv );
  // Assume arguments are sanitized
  
  // get the input file
  auto mesh_filename = 
    args.count("m") ? args.at("m") : std::string();
  
  // need to put the filename into a statically sized character array
  auto mesh_basename = utils::basename( mesh_filename );

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // get the client handle 
  auto mesh = flecsi_get_client_handle(mesh_t, meshes, mesh0);
 
  // cout << mesh;
  
  // output the mesh
  auto prefix = utils::remove_extension( mesh_basename );
  char_array_t prefix_char;
  strcpy( prefix_char.data(), prefix.c_str() );
  auto f1 = flecsi_execute_task(output, single, mesh, prefix_char);
  f1.wait();

  return 0;

#if 0

  //===========================================================================
  // Some typedefs
  //===========================================================================

  using size_t = typename mesh_t::size_t;
  using real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t; 

  // get machine zero
  constexpr auto epsilon = std::numeric_limits<real_t>::epsilon();
  const auto machine_zero = std::sqrt(epsilon);

  // the maximum number of retries
  constexpr int max_retries = 5;

  //===========================================================================
  // Field Creation
  //===========================================================================

  // start the timer
  auto tstart = utils::get_wall_time();



  // set these variables as persistent for plotting
  flecsi_get_accessor(mesh, hydro,  density,   real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, pressure,   real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, velocity, vector_t, dense, 0).attributes().set(persistent);

  flecsi_get_accessor(mesh, hydro, internal_energy, real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro,     temperature, real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro,     sound_speed, real_t, dense, 0).attributes().set(persistent);
#endif

  //===========================================================================
  // Access what we need
  //===========================================================================
  
  auto d  = flecsi_get_handle(mesh, hydro,  density,   real_t, dense, 0);

#if 0
  auto d0 = flecsi_get_handle(mesh, hydro,  density,   real_t, dense, 1);
  auto v  = flecsi_get_handle(mesh, hydro, velocity, vector_t, dense, 0);
  auto v0 = flecsi_get_handle(mesh, hydro, velocity, vector_t, dense, 1);
  auto e  = flecsi_get_handle(mesh, hydro, internal_energy, real_t, dense, 0);
  auto e0 = flecsi_get_handle(mesh, hydro, internal_energy, real_t, dense, 1);

  auto p  = flecsi_get_handle(mesh, hydro,        pressure,   real_t, dense, 0);
  auto T  = flecsi_get_handle(mesh, hydro,     temperature, real_t, dense, 0);
  auto a  = flecsi_get_handle(mesh, hydro,     sound_speed, real_t, dense, 0);

  auto F = flecsi_get_handle(mesh, hydro, flux, flux_data_t, dense, 0);
#endif

  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  flecsi_execute_task( 
    initial_conditions, 
    single, 
    mesh, 
    inputs_t::ics,
    d
  );

#if 0
 
  #ifdef HAVE_CATALYST
    auto insitu = io::catalyst::adaptor_t(catalyst_scripts);
    std::cout << "Catalyst on!" << std::endl;
  #endif

  // Update the EOS
  flecsi_execute_task( 
    update_state_from_pressure_task, loc, single, mesh, inputs_t::eos.get() 
  );

  //===========================================================================
  // Pre-processing
  //===========================================================================


  // now output the solution
  if (inputs_t::output_freq > 0)
    output(mesh, inputs_t::prefix, inputs_t::postfix, 1);

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  // get the current time
  auto soln_time = mesh.time();
  auto time_cnt  = mesh.time_step_counter();

  // get an accessor for the time step
  auto time_step = flecsi_get_accessor( mesh, hydro, time_step, real_t, global, 0 );   

  // a counter for this session
  size_t num_steps = 0; 

  for ( size_t num_retries = 0;
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    // store the initial solution, only if this isnt a retry
    if (num_retries == 0)
      flecsi_execute_task( save_solution_task, loc, single, mesh );

    // compute the time step
    flecsi_execute_task( evaluate_time_step_task, loc, single, mesh );
 
    // access the computed time step and make sure its not too large
    *time_step = std::min( *time_step, inputs_t::final_time - soln_time );       

    //-------------------------------------------------------------------------
    // try a timestep

    // compute the fluxes
    flecsi_execute_task( evaluate_fluxes_task, loc, single, mesh );

    // reset the time stepping mode
    auto mode = mode_t::normal;

    // wrap the stage update in a do loop, so it can be re-executed
    do {

      // output the time step
      cout << std::string(80, '=') << endl;
      auto ss = cout.precision();
      cout.setf( std::ios::scientific );
      cout.precision(6);
      cout << "|  " << "Step:" << std::setw(10) << time_cnt+1
           << "  |  Time:" << std::setw(17) << soln_time + (*time_step)
           << "  |  Step Size:" << std::setw(17) << *time_step
           << "  |" << std::endl;
      cout.unsetf( std::ios::scientific );
      cout.precision(ss);

      // Loop over each cell, scattering the fluxes to the cell
      auto err = 
        flecsi_execute_task( 
          apply_update_task, loc, single, mesh, machine_zero, true 
        );
      auto update_flag = err.get();


      // dump the current errored solution to a file
      if ( update_flag != solution_error_t::ok && inputs_t::output_freq > 0)
        output(mesh, inputs_t::prefix+"-error", inputs_t::postfix, 1);

      // if we got an unphysical solution, half the time step and try again
      if ( update_flag == solution_error_t::unphysical ) {
        // Print a message that we are retrying
        std::cout << "Unphysical solution detected, halfing timestep..." << std::endl;
        // half time step
        *time_step *= 0.5;
        // indicate that we are retrying
        mode = mode_t::retry;
      }

      // if there was variance, retry the solution again
      else if (update_flag == solution_error_t::variance) {
        // Print a message that we are retrying
        std::cout << "Variance in solution detected, retrying..." << std::endl;
        // roll the counter back
        --num_steps;
        // set the flag to restart
        mode = mode_t::restart;
      }

      // if there is no error, just break
      else {
        mode = mode_t::normal;
        break;
      }


      // if we are retrying or restarting, restore the original solution
      if (mode==mode_t::retry || mode==mode_t::restart) {
        // restore the initial solution
        flecsi_execute_task( restore_solution_task, loc, single, mesh );
        // don't retry forever
        if ( ++num_retries > max_retries ) {
          // Print a message we are exiting
          std::cout << "Too many retries, exiting..." << std::endl;
          // flag that we want to quit
          mode = mode_t::quit;
        }
      }


    } while ( mode==mode_t::retry );

    // end timestep
    //-------------------------------------------------------------------------

    // if a restart is detected, restart the whole iteration loop
    if (mode==mode_t::restart) continue;

    // Update derived solution quantities
    flecsi_execute_task( 
      update_state_from_energy_task, loc, single, mesh, inputs_t::eos.get() 
    );

    // now we can quit after the solution has been reset to the previous step's
    if (mode==mode_t::quit) break;

    #ifdef HAVE_CATALYST
    if (!catalyst_scripts.empty()) {
      auto vtk_grid = mesh::to_vtk( mesh );
      insitu.process( 
        vtk_grid, soln_time, num_steps, (num_steps==inputs_t::max_steps-1)
      );
    }
    #endif

    // update time
    soln_time = mesh.increment_time( *time_step );
    time_cnt = mesh.increment_time_step_counter();

    // now output the solution
    output(mesh, inputs_t::prefix, inputs_t::postfix, inputs_t::output_freq);

    // reset the number of retrys if we eventually made it through a time step
    num_retries  = 0;

  }

  //===========================================================================
  // Post-process
  //===========================================================================
    
  // now output the solution
  if ( (inputs_t::output_freq > 0) && (time_cnt % inputs_t::output_freq != 0) )
    output(mesh, inputs_t::prefix, inputs_t::postfix, 1);

  cout << "Final solution time is " 
       << std::scientific << std::setprecision(2) << soln_time
       << " after " << num_steps << " steps." << std::endl;

  
  auto tdelta = utils::get_wall_time() - tstart;
  std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
            << tdelta << "s." << std::endl;


  // now output the checksums
  mesh::checksum(mesh);

#endif

  // success if you reached here
  return 0;

}

} // namespace
} // namespace
