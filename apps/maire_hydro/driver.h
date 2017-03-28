/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro incdludes
#include "types.h"
#include "../common/exceptions.h"
#include "../common/parse_arguments.h"

// user includes
#include <flecsale/eos/ideal_gas.h>
#include <flecsale/mesh/mesh_utils.h>
#include <flecsale/utils/time_utils.h>

// system includes
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

namespace apps {
namespace hydro {

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
template< typename inputs_t >
int driver(int argc, char** argv) 
{


  // set exceptions 
  enable_exceptions();

  //===========================================================================
  // Parse arguments
  //===========================================================================

  // the usage stagement
  auto print_usage = [&argv]() {
    std::cout << "Usage: " << argv[0] 
              << " [--file INPUT_FILE]"
              << " [--help]"
              << std::endl << std::endl;
    std::cout << "\t--file INPUT_FILE:\t Override the input file "
              << "with INPUT_FILE." << std::endl;
    std::cout << "\t--help:\t Print a help message." << std::endl;
  };

  // Define the options
  struct option long_options[] =
    {
      {"help",       no_argument, 0, 'h'},
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };
  const char * short_options = "hf:";

  // parse the arguments
  auto args = parse_arguments(argc, argv, long_options, short_options);

  // process the simple ones
  if ( args.count("h") ) {
    print_usage();
    return 0;
  } 
  else if ( args.count("?") ) {
    print_usage();
    return 1;
  }
 
  // get the input file
  auto input_file_name = 
    args.count("f") ? args.at("f") : std::string();

  // override any inputs if need be
  if ( !input_file_name.empty() ) {
    std::cout << "Using input file \"" << input_file_name << "\"." 
              << std::endl;
    inputs_t::load( input_file_name );
  }

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // make the mesh
  auto mesh = inputs_t::make_mesh( /* solution time */ 0.0 );

  // this is the mesh object
  mesh.is_valid();
  
  cout << mesh;

  //===========================================================================
  // Some typedefs
  //===========================================================================

  using mesh_t = decltype(mesh);
  using size_t = typename mesh_t::size_t;
  using real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t; 

  // get machine zero
  constexpr auto epsilon = std::numeric_limits<real_t>::epsilon();
  const auto machine_zero = std::sqrt(epsilon);

  // the maximum number of retries
  constexpr int max_retries = 5;

  // the time stepping stage coefficients
  std::vector<real_t> stages = {0.5, 1.0};

  //===========================================================================
  // Field Creation
  //===========================================================================

  // start the timer
  auto tstart = utils::get_wall_time();

  // type aliases
  using matrix_t = matrix_t< mesh_t::num_dimensions >;
  using flux_data_t = flux_data_t< mesh_t::num_dimensions >;

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  flecsi_register_data(mesh, hydro, cell_volume,     real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro, cell_mass,       real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro, cell_pressure,   real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro, cell_velocity, vector_t, dense, 2, cells);

  flecsi_register_data(mesh, hydro, cell_density,         real_t, dense, 2, cells);
  flecsi_register_data(mesh, hydro, cell_internal_energy, real_t, dense, 2, cells);
  flecsi_register_data(mesh, hydro, cell_temperature,     real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro, cell_sound_speed,     real_t, dense, 1, cells);

  // node state
  flecsi_register_data(mesh, hydro, node_coordinates, vector_t, dense, 1, vertices);
  flecsi_register_data(mesh, hydro, node_velocity, vector_t, dense, 1, vertices);

  // solver state
  flecsi_register_data(mesh, hydro, cell_residual, flux_data_t, dense, 1, cells);

  flecsi_register_data(mesh, hydro, corner_normal, vector_t, dense, 1, corners);
  flecsi_register_data(mesh, hydro, corner_force, vector_t, dense, 1, corners);
  
  // register the time step and set a cfl
  flecsi_register_data( mesh, hydro, time_step, real_t, global, 1 );
  flecsi_register_data( mesh, hydro, cfl, time_constants_t, global, 1 );
  
  *flecsi_get_accessor( mesh, hydro, time_step, real_t, global, 0 ) = inputs_t::initial_time_step;
  *flecsi_get_accessor( mesh, hydro, cfl, time_constants_t, global, 0 ) = inputs_t::CFL;

  // Register the total energy
  flecsi_register_data( mesh, hydro, sum_total_energy, real_t, global, 1 );

  // set the persistent variables, i.e. the ones that will be plotted
  flecsi_get_accessor(mesh, hydro, cell_mass,       real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, cell_pressure,   real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, cell_velocity, vector_t, dense, 0).attributes().set(persistent);

  flecsi_get_accessor(mesh, hydro, cell_density,         real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, cell_internal_energy, real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, cell_temperature,     real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, cell_sound_speed,     real_t, dense, 0).attributes().set(persistent);

  flecsi_get_accessor(mesh, hydro, node_velocity, vector_t, dense, 0).attributes().set(persistent);


  //===========================================================================
  // Boundary Conditions
  //===========================================================================
  
  // get the current time
  auto soln_time = mesh.time();
  auto time_cnt  = mesh.time_step_counter();

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // install each boundary
  for ( const auto & bc_pair : inputs_t::bcs )
  {
    auto bc_type = bc_pair.first.get();
    auto bc_function = bc_pair.second; 
    auto bc_key = mesh.install_boundary( 
      [=](auto f) 
      { 
        if ( f->is_boundary() ) {
          const auto & fx = f->midpoint();
          return ( bc_function(fx, soln_time) );
        }
        return false;
      }
    );
    boundaries.emplace( bc_key, bc_type );
  }


  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  flecsi_execute_task( initial_conditions_task, loc, single, mesh, inputs_t::ics );
  

  // Update the EOS
  flecsi_execute_task( 
    update_state_from_pressure_task, loc, single, mesh, inputs_t::eos.get()
  );


  //===========================================================================
  // Pre-processing
  //===========================================================================

  // now output the solution
  if ( inputs_t::output_freq > 0 )
    output(mesh, inputs_t::prefix, inputs_t::postfix, 1);
  

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  // get the time step accessor
  auto time_step = flecsi_get_accessor( mesh, hydro, time_step, real_t, global, 0 );   

  // a counter for this session
  size_t num_steps = 0; 

  for (
    size_t num_retries = 0;
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    //--------------------------------------------------------------------------
    // Begin Time step
    //--------------------------------------------------------------------------

    // Save solution at n=0
    if (num_retries == 0) {
      flecsi_execute_task( save_coordinates_task, loc, single, mesh );
      flecsi_execute_task( save_solution_task, loc, single, mesh );
    }

    // keep the old time step
    real_t time_step_old = time_step;

    //--------------------------------------------------------------------------
    // Predictor step : Evaluate Forces at n=0
    //--------------------------------------------------------------------------

    // estimate the nodal velocity at n=0
    flecsi_execute_task( estimate_nodal_state_task, loc, single, mesh );

    // compute the nodal velocity at n=0
    flecsi_execute_task( 
      evaluate_nodal_state_task, loc, single, mesh, boundaries
    );

    // compute the fluxes
    flecsi_execute_task( evaluate_residual_task, loc, single, mesh );

    //--------------------------------------------------------------------------
    // Time step evaluation
    //--------------------------------------------------------------------------

    // compute the time step
    std::string limit_string;
    flecsi_execute_task( evaluate_time_step_task, loc, single, mesh, limit_string );
    
    // access the computed time step and make sure its not too large
    *time_step = std::min( *time_step, inputs_t::final_time - soln_time );       

    cout << std::string(60, '=') << endl;
    auto ss = cout.precision();
    cout.setf( std::ios::scientific );
    cout.precision(6);
    cout << "| " << std::setw(8) << "Step:"
         << " | " << std::setw(13) << "Time:"
         << " | " << std::setw(13) << "Step Size:"
         << " | " << std::setw(13) << "Limit:"
         << " |" << std::endl;
    cout << "| " << std::setw(8) << time_cnt+1
         << " | " << std::setw(13) << soln_time + (*time_step)
         << " | " << std::setw(13) << *time_step
         << " | " << std::setw(13) << limit_string
         << " |" << std::endl;
    cout.unsetf( std::ios::scientific );
    cout.precision(ss);

    //--------------------------------------------------------------------------
    // Multistage 
    //--------------------------------------------------------------------------
    
    // a stage counter
    int istage = 0;

    // reset the time stepping mode
    auto mode = mode_t::normal;

    do {

      //------------------------------------------------------------------------
      // Move to n^stage

      // move the mesh to n+1/2
      flecsi_execute_task( move_mesh_task, loc, single, mesh, stages[istage] );

      // update solution to n+1/2
      auto err = flecsi_execute_task( 
        apply_update_task, loc, single, mesh, stages[istage], machine_zero, (istage==0)
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
        // reset time step
        *time_step = time_step_old;
        // roll the counter back
        --num_steps;
        // set the flag to restart
        mode = mode_t::restart;
      }

      // otherewise, set the mode to none
      else {
        mode = mode_t::normal;
      }

      // if we are retrying or restarting, restore the original solution
      if (mode == mode_t::restart || mode == mode_t::retry) {
        // restore the initial solution
        flecsi_execute_task( restore_coordinates_task, loc, single, mesh );
        flecsi_execute_task( restore_solution_task, loc, single, mesh );
        mesh.update_geometry();
        // don't retry forever
        if ( ++num_retries > max_retries ) {
          // Print a message we are exiting
          std::cout << "Too many retries, exiting..." << std::endl;
          // flag that we want to quit
          mode = mode_t::quit;
        }
      }

      // Update derived solution quantities
      flecsi_execute_task( 
        update_state_from_energy_task, loc, single, mesh, inputs_t::eos.get() 
      );

      // compute the current nodal velocity
      flecsi_execute_task( 
        evaluate_nodal_state_task, loc, single, mesh, boundaries
      );

      // if we are retrying, then restart the loop since all the state has been 
      // reset
      if (mode == mode_t::retry) {
        istage = 0;
        continue;
      }

      // now we can quit once all the state has been reset
      else if ( 
        mode == mode_t::quit || 
        mode == mode_t::restart || 
        ++istage==stages.size()
      ) 
        break;

      //------------------------------------------------------------------------
      // Corrector : Evaluate Forces at n^stage

      // compute the fluxes
      flecsi_execute_task( evaluate_residual_task, loc, single, mesh );

      //------------------------------------------------------------------------
      // Move to n+1

      // restore the solution to n=0
      flecsi_execute_task( restore_coordinates_task, loc, single, mesh );
      flecsi_execute_task( restore_solution_task, loc, single, mesh );

    } while(true); // do

    //--------------------------------------------------------------------------
    // End Time step
    //--------------------------------------------------------------------------

    // now we can quit once all the state has been reset
    if (mode == mode_t::quit) break;
    // now we can restart once all the state has been reset
    else if (mode == mode_t::restart) continue;

    // update time
    soln_time = mesh.increment_time( *time_step );
    time_cnt = mesh.increment_time_step_counter();
  
    // now output the solution
    output(
      mesh, inputs_t::prefix, inputs_t::postfix, inputs_t::output_freq
    );

    // if we got through a whole cycle, reset the retry counter
    num_retries = 0;

  }


  //===========================================================================
  // Post-process
  //===========================================================================
    
  // now output the solution
  if ( (inputs_t::output_freq > 0) && (time_cnt % inputs_t::output_freq != 0) )
    output(mesh, inputs_t::prefix, inputs_t::postfix, 1);

  cout << "Final solution time is " 
       << std::scientific << std::setprecision(6) << soln_time
       << " after " << num_steps << " steps." << std::endl;

  
  auto tdelta = utils::get_wall_time() - tstart;
  std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
            << tdelta << "s." << std::endl;
  
  // now output the checksums
  mesh::checksum(mesh);

  // success
  return 0;

}

} // namespace
} // namespace
