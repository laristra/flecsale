/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "arguments.h"
#include "types.h"

// user includes
//#include <flecsale/mesh/mesh_utils.h>
#include <flecsale/utils/time_utils.h>

// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

namespace apps {
namespace hydro {
  
// register the data client
flecsi_register_data_client(mesh_t, meshes, mesh0);

// create some field data.  Fields are registered as struct of arrays.
// this allows us to access the data in different patterns.

// The cell state
flecsi_register_field(
  mesh_t, 
  hydro,  
  cell_volume,
  mesh_t::real_t, 
  dense, 
  1, 
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro,  
  cell_mass,
  mesh_t::real_t, 
  dense, 
  1, 
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro, 
  cell_pressure,
  mesh_t::real_t, 
  dense, 
  1, 
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro, 
  cell_velocity,
  mesh_t::vector_t,
  dense,
  2,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro,  
  cell_density,   
  mesh_t::real_t, 
  dense, 
  2, 
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro,
  cell_internal_energy,
  mesh_t::real_t,
  dense,
  2,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,
  hydro,
  cell_temperature,
  mesh_t::real_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,
  hydro,
  cell_sound_speed,
  mesh_t::real_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

// node state
flecsi_register_field(
  mesh_t,
  hydro,
  node_coordinates,
  mesh_t::vector_t,
  dense,
  1,
  mesh_t::index_spaces_t::vertices
);

flecsi_register_field(
  mesh_t,
  hydro,
  node_velocity,
  mesh_t::vector_t,
  dense,
  1,
  mesh_t::index_spaces_t::vertices
);

// solver state
flecsi_register_field(
  mesh_t,
  hydro,
  cell_residual,
  flux_data_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,
  hydro,
  corner_normal,
  mesh_t::vector_t,
  dense,
  1,
  mesh_t::index_spaces_t::corners
);

flecsi_register_field(
  mesh_t,
  hydro,
  corner_force,
  mesh_t::vector_t,
  dense,
  1,
  mesh_t::index_spaces_t::corners
);
  

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int driver(int argc, char** argv) 
{


  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

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

  // figure out an output prefix, this will be used later
  auto prefix = utils::remove_extension( mesh_basename );

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // start a clock
  auto tstart = utils::get_wall_time();

  // get the client handle 
  auto mesh = flecsi_get_client_handle(mesh_t, meshes, mesh0);
 
  // cout << mesh;
  
  //===========================================================================
  // Some typedefs
  //===========================================================================

  using size_t = typename mesh_t::size_t;
  using real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t; 

  // get machine zero
  constexpr auto epsilon = std::numeric_limits<real_t>::epsilon();
  const auto machine_zero = std::sqrt(epsilon);

  // the time stepping stage coefficients
  std::vector<real_t> stages = {0.5, 1.0};
  
  // the solution time starts at zero
  real_t soln_time{0};  
  size_t time_cnt{0}; 

  //===========================================================================
  // Access what we need
  //===========================================================================

  // cell data
  auto Vc = flecsi_get_handle(mesh, hydro, cell_volume,     real_t, dense, 0);
  auto Mc = flecsi_get_handle(mesh, hydro, cell_mass,       real_t, dense, 0);
  auto pc = flecsi_get_handle(mesh, hydro, cell_pressure,   real_t, dense, 0);
  auto uc = flecsi_get_handle(mesh, hydro, cell_velocity, vector_t, dense, 0);

  auto dc = flecsi_get_handle(mesh, hydro, cell_density,         real_t, dense, 0);
  auto ec = flecsi_get_handle(mesh, hydro, cell_internal_energy, real_t, dense, 0);
  auto Tc = flecsi_get_handle(mesh, hydro, cell_temperature,     real_t, dense, 0);
  auto ac = flecsi_get_handle(mesh, hydro, cell_sound_speed,     real_t, dense, 0);

  // node state
  auto xn = flecsi_get_handle(mesh, hydro, node_coordinates, vector_t, dense, 0);
  auto un = flecsi_get_handle(mesh, hydro, node_velocity, vector_t, dense, 0);

  // solver state
  auto dUdt = flecsi_get_handle(mesh, hydro, cell_residual, flux_data_t, dense, 0);
  auto npc = flecsi_get_handle(mesh, hydro, corner_normal, vector_t, dense, 0);
  auto Fpc = flecsi_get_handle(mesh, hydro, corner_force, vector_t, dense, 0);
  

  //===========================================================================
  // Boundary Conditions
  //===========================================================================
  

#if 0
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
    inputs_t::eos,
    soln_time,
    Vc, Mc, uc, pc, dc, ec, Tc, ac
  );


  //===========================================================================
  // Pre-processing
  //===========================================================================


  // now output the solution
  auto has_output = (inputs_t::output_freq > 0);
  if (has_output) {
    auto name = prefix + "_" + apps::common::zero_padded( 0 ) + ".exo";
    auto name_char = utils::to_trivial_string( name );
    auto f =
      flecsi_execute_task(
        output, single, mesh, name_char, dc, uc, ec, pc, Tc, ac
      );
    f.wait();
  }



  // dump connectivity
  {
    auto name = utils::to_trivial_string( prefix+".txt" );
    auto f = flecsi_execute_task(print, single, mesh, name);
    f.wait();
  }

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  // Get the legion runtime and context.  In the furture, legion will not be 
  // exposed to the user
  auto legion_runtime = Legion::Runtime::get_runtime();
  auto legion_context = Legion::Runtime::get_context();

  auto & min_reduction = 
    flecsi::execution::context_t::instance().min_reduction();

  for (
    size_t num_steps = 0;
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    //--------------------------------------------------------------------------
    // Predictor step : Evaluate Forces at n=0
    //--------------------------------------------------------------------------

    // estimate the nodal velocity at n=0
    flecsi_execute_task( estimate_nodal_state, single, mesh, uc, un );

    // compute the nodal velocity at n=0
    flecsi_execute_task( 
      evaluate_nodal_state,
      single,
      mesh,
      soln_time,
      Vc, Mc, uc, pc, dc, ec, Tc, ac,
      un, npc, Fpc
    );

#if 0
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

#endif

    // now output the solution
    if ( has_output && 
        (time_cnt % inputs_t::output_freq == 0 || 
         num_steps==inputs_t::max_steps-1 ||
         std::abs(soln_time-inputs_t::final_time) < epsilon
        )  
      ) 
    {
      auto name = prefix + "_" + apps::common::zero_padded( time_cnt ) + ".exo";
      auto name_char = utils::to_trivial_string( name );
      auto f = flecsi_execute_task(
        output, single, mesh, name_char, dc, uc, ec, pc, Tc, ac
      );
      f.wait();
    }

  }


  //===========================================================================
  // Post-process
  //===========================================================================

  auto tdelta = utils::get_wall_time() - tstart;

  if ( rank == 0 ) {

    cout << "Final solution time is " 
         << std::scientific << std::setprecision(2) << soln_time
         << " after " << time_cnt << " steps." << std::endl;

    
    std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
              << tdelta << "s." << std::endl;

  }


  // success if you reached here
  return 0;

}

} // namespace
} // namespace
