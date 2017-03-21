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
#include "types.h"
#include "../common/exceptions.h"
#include "../common/parse_arguments.h"

// user includes
#include <flecsale/mesh/mesh_utils.h>
#include <flecsale/utils/time_utils.h>

// system includes
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#ifdef USE_CATALYST
  #include "flecsaleAdaptor.h"
  #include "flecsale_unstructuredGrid.h"
#endif

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
  flecsi_register_data(mesh, hydro,  density,   real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro, pressure,   real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro, velocity, vector_t, dense, 1, cells);

  flecsi_register_data(mesh, hydro, internal_energy, real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro,     temperature, real_t, dense, 1, cells);
  flecsi_register_data(mesh, hydro,     sound_speed, real_t, dense, 1, cells);

  // set these variables as persistent for plotting
  flecsi_get_accessor(mesh, hydro,  density,   real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, pressure,   real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro, velocity, vector_t, dense, 0).attributes().set(persistent);

  flecsi_get_accessor(mesh, hydro, internal_energy, real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro,     temperature, real_t, dense, 0).attributes().set(persistent);
  flecsi_get_accessor(mesh, hydro,     sound_speed, real_t, dense, 0).attributes().set(persistent);

  // compute the fluxes.  here I am regestering a struct as the stored data
  // type since I will only ever be accesissing all the data at once.
  flecsi_register_data(mesh, hydro, flux, flux_data_t, dense, 1, faces);

  // register the time step and set a cfl
  flecsi_register_data( mesh, hydro, time_step, real_t, global, 1 );
  flecsi_register_data( mesh, hydro, cfl, real_t, global, 1 );
  *flecsi_get_accessor( mesh, hydro, cfl, real_t, global, 0) = inputs_t::CFL;  


  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  flecsi_execute_task( initial_conditions_task, loc, single, mesh, inputs_t::ics );
  
  #ifdef USE_CATALYST
    FlecsaleAdaptor::Init(argc, argv);
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

  for ( ; 
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    // compute the time step
    flecsi_execute_task( evaluate_time_step_task, loc, single, mesh );
 
    // access the computed time step and make sure its not too large
    *time_step = std::min( *time_step, inputs_t::final_time - soln_time );       

    cout << "step =  " << std::setw(4) << time_cnt+1
         << std::setprecision(2)
         << ", time = " << std::scientific << soln_time + (*time_step)
         << ", dt = " << std::scientific << *time_step
         << std::endl;

    // compute the fluxes
    flecsi_execute_task( evaluate_fluxes_task, loc, single, mesh );
    
    // Loop over each cell, scattering the fluxes to the cell
    flecsi_execute_task( apply_update_task, loc, single, mesh );

    // Update derived solution quantities
    flecsi_execute_task( 
      update_state_from_energy_task, loc, single, mesh, inputs_t::eos.get() 
    );


    #ifdef USE_CATALYST
      vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
      populate(mesh, grid, num_steps);

      FlecsaleAdaptor::CoProcess( grid, soln_time, num_steps, (num_steps==inputs_t::max_steps-1) );

      // update time
      soln_time = mesh.increment_time( *time_step );
      time_cnt = mesh.increment_time_step_counter();

      //std::cout << "fed to catalyst...\n" << std::endl;
    #else
      // update time
      soln_time = mesh.increment_time( *time_step );
      time_cnt = mesh.increment_time_step_counter();

      // now output the solution
      output(mesh, inputs_t::prefix, inputs_t::postfix, inputs_t::output_freq);
    #endif

   

  }

  //===========================================================================
  // Post-process
  //===========================================================================
    
  #ifdef USE_CATALYST
    FlecsaleAdaptor::Finalize();
  #endif   
    
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


  // success if you reached here
  return 0;

}

} // namespace
} // namespace
