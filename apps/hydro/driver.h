/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This is the main driver for the hydro solver.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "tasks.h"
#include "types.h"
#include "../common/exceptions.h"
#include "../common/parse_arguments.h"

// user includes
#include <ale/utils/time_utils.h>

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
  auto mesh = inputs_t::make_mesh();

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

  //===========================================================================
  // Field Creation
  //===========================================================================

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  // start the timer
  auto tstart = utils::get_wall_time();


  // type aliases
  using eqns_t = eqns_t<mesh_t::num_dimensions>;
  using flux_data_t = flux_data_t<mesh_t::num_dimensions>;

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  register_data(mesh, hydro,  density,   real_t, dense, 1, cells);
  register_data(mesh, hydro, pressure,   real_t, dense, 1, cells);
  register_data(mesh, hydro, velocity, vector_t, dense, 1, cells);

  register_data(mesh, hydro, internal_energy, real_t, dense, 1, cells);
  register_data(mesh, hydro,     temperature, real_t, dense, 1, cells);
  register_data(mesh, hydro,     sound_speed, real_t, dense, 1, cells);

  // set these variables as persistent for plotting
  get_accessor(mesh, hydro,  density,   real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, pressure,   real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, velocity, vector_t, dense, 0).attributes().set(persistent);

  get_accessor(mesh, hydro, internal_energy, real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro,     temperature, real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro,     sound_speed, real_t, dense, 0).attributes().set(persistent);

  // compute the fluxes.  here I am regestering a struct as the stored data
  // type since I will only ever be accesissing all the data at once.
  register_data(mesh, hydro, flux, flux_data_t, dense, 1, faces);

  // register the time step and set a cfl
  register_data( mesh, hydro, time_step, real_t, global, 1 );
  register_data( mesh, hydro, cfl, real_t, global, 1 );
  *get_accessor( mesh, hydro, cfl, real_t, global, 0) = inputs_t::CFL;  

  // register state a global eos
  register_data( mesh, hydro, eos, eos_t, global, 1 );
  *get_accessor( mesh, hydro, eos, eos_t, global, 0 ) = eos;

  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  initial_conditions( mesh, inputs_t::ics );
  

  // Update the EOS
  update_state_from_pressure( mesh );

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
  auto time_step = get_accessor( mesh, hydro, time_step, real_t, global, 0 );   

  // a counter for this session
  size_t num_steps = 0; 

  for ( ; 
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    // compute the time step
    evaluate_time_step<eqns_t>( mesh );
    
    // access the computed time step and make sure its not too large
    *time_step = std::min( *time_step, inputs_t::final_time - soln_time );       

    cout << "step =  " << std::setw(4) << time_cnt+1
         << std::setprecision(2)
         << ", time = " << std::scientific << soln_time + (*time_step)
         << ", dt = " << std::scientific << *time_step
         << std::endl;

    // compute the fluxes
    evaluate_fluxes( mesh );
    
    // Loop over each cell, scattering the fluxes to the cell
    apply_update( mesh );

    // Update derived solution quantities
    update_state_from_energy( mesh );

    // update time
    soln_time = mesh.increment_time( *time_step );
    time_cnt = mesh.increment_time_step_counter();

    // now output the solution
    output(mesh, inputs_t::prefix, inputs_t::postfix, inputs_t::output_freq);

  }

  //===========================================================================
  // Post-process
  //===========================================================================
    
  // now output the solution
  if ( (inputs_t::output_freq > 0) && (time_cnt%inputs_t::output_freq != 0) )
    output(mesh, inputs_t::prefix, inputs_t::postfix, 1);

  cout << "Final solution time is " 
       << std::scientific << std::setprecision(2) << soln_time
       << " after " << num_steps << " steps." << std::endl;

  
  auto tdelta = utils::get_wall_time() - tstart;
  std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
            << tdelta << "s." << std::endl;

  // success if you reached here
  return 0;

}

} // namespace
} // namespace
