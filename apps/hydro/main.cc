/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro incdludes
#include "tasks.h"
#include "types.h"
#include "../common/exceptions.h"

// user includes
#include <ale/mesh/factory.h>
#include <ale/utils/python_utils.h>
#include <ale/utils/time_utils.h>

// system includes
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

// everything is in the hydro namespace
using namespace apps::hydro;

// right now cases are hard coded
//#define SODX_2D
//#define SODX_3D
//#define SHOCK_BOX_2D
#define SHOCK_BOX_3D


///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

  // set exceptions 
  enable_exceptions();

  //===========================================================================
  // Parse arguments
  //===========================================================================

  // argument defaults
  bool is_verbose = false;
  std::string input_file_name = "infile.py";

  // the usage stagement
  auto print_usage = [&argv]() {
    std::cout << "Usage: " << argv[0] 
              << " [--verbose]"
              << " [--help]"
              << " [--file INPUT_FILE]"
              << std::endl << std::endl;
    std::cout << "\t--file INPUT_FILE:\t Override the input file "
              << "with INPUT_FILE." << std::endl;
    std::cout << "\t--help:\t Print the usage." << std::endl;
    std::cout << "\t--verbose:\t Turn on verbosity." << std::endl;
  };

  while (1) {
    
    // Define the options
    static struct option long_options[] =
      {
        {"help",            no_argument, 0, 'h'},
        {"verbose",         no_argument, 0, 'v'},
        {"file",      required_argument, 0, 'f'},
        {0, 0, 0, 0}
      };
      // getopt_long stores the option index here.
      int option_index = 0;

      auto c = getopt_long (argc, argv, "hvf:", long_options, &option_index);

      // Detect the end of the options.
      if (c == -1) break;

      switch (c) {
        case 'f':
          short_warning ("Using input file \"" << optarg << "\"." );
          input_file_name = optarg;
          break;

        case 'v':
          short_warning ("Setting verbosity on." );
          is_verbose = true;
          break;

        case 'h':
          print_usage();
          return 0;

        case '?':
          // getopt_long already printed an error message.
          print_usage();
          return 1;

        default:
          abort ();
        }
    }

  //===========================================================================
  // Case Setup
  //===========================================================================

  // setup the python interpreter
  utils::python_set_program_name( argv[0] );
  auto py_module = utils::python_import( input_file_name.c_str() );

  utils::python_initialize();

  // shut down python
  utils::python_finalize();

#if 0


  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // this is the mesh object
  mesh.is_valid();
  
  cout << mesh;

  //===========================================================================
  // Field Creation
  //===========================================================================

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
  *get_accessor( mesh, hydro, cfl, real_t, global, 0) = CFL;  

  // register state a global eos
  register_data( mesh, hydro, eos, eos_t, global, 1 );
  *get_accessor( mesh, hydro, eos, eos_t, global, 0 ) = eos;

  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  apps::hydro::initial_conditions( mesh, ics );
  

  // Update the EOS
  apps::hydro::update_state_from_pressure( mesh );

  //===========================================================================
  // Pre-processing
  //===========================================================================


  // now output the solution
  if (output_freq > 0)
    apps::hydro::output(mesh, prefix, postfix, 1);

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

  for ( ; (num_steps < max_steps && soln_time < final_time); ++num_steps ) 
  {   

    // compute the time step
    apps::hydro::evaluate_time_step<eqns_t>( mesh );
    
    // access the computed time step and make sure its not too large
    *time_step = std::min( *time_step, final_time - soln_time );       

    cout << "step =  " << std::setw(4) << time_cnt+1
         << std::setprecision(2)
         << ", time = " << std::scientific << soln_time + (*time_step)
         << ", dt = " << std::scientific << *time_step
         << std::endl;

    // compute the fluxes
    apps::hydro::evaluate_fluxes( mesh );
    
    // Loop over each cell, scattering the fluxes to the cell
    apps::hydro::apply_update( mesh );

    // Update derived solution quantities
    apps::hydro::update_state_from_energy( mesh );

    // update time
    soln_time = mesh.increment_time( *time_step );
    time_cnt = mesh.increment_time_step_counter();

    // now output the solution
    apps::hydro::output(mesh, prefix, postfix, output_freq);

  }

  //===========================================================================
  // Post-process
  //===========================================================================
    
  // now output the solution
  if ( (output_freq > 0) && (time_cnt%output_freq != 0) )
    apps::hydro::output(mesh, prefix, postfix, 1);

  cout << "Final solution time is " 
       << std::scientific << std::setprecision(2) << soln_time
       << " after " << num_steps << " steps." << std::endl;

  
  auto tdelta = utils::get_wall_time() - tstart;
  std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
            << tdelta << "s." << std::endl;

#endif 

  // success
  return 0;

}

