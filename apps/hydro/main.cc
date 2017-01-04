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

// user includes
#include <ale/mesh/factory.h>
#include <ale/utils/time_utils.h>

// system includes
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#ifdef _GNU_SOURCE
#  include <fenv.h>
#endif

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

  // enable exceptions
#if defined(_GNU_SOURCE) && !defined(NDEBUG)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  //===========================================================================
  // X-direction SOD Inputs
  //===========================================================================

#ifdef SODX_2D

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sodx_2d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 1;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr real_t CFL = 1.0;
  constexpr real_t final_time = 0.2;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[0] < 0.0 ) {
        d = 1.0;
        p = 1.0;
      }
      else {
        d = 0.125;
        p = 0.1;
      }    
      return std::make_tuple( d, v, p );
    };

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  //===========================================================================
  // X-direction SOD Inputs
  //===========================================================================
#elif defined(SODX_3D)

  using mesh_t = mesh_3d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sodx_3d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 1;
  constexpr size_t num_cells_z = 1;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  constexpr real_t length_z = 0.1;
  
  // the CFL and final solution time
  constexpr real_t CFL = 1.0;
  constexpr real_t final_time = 0.2;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[0] < 0.0 ) {
        d = 1.0;
        p = 1.0;
      }
      else {
        d = 0.125;
        p = 0.1;
      }    
      return std::make_tuple( d, v, p );
    };

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( 
    num_cells_x, num_cells_y, num_cells_z, length_x, length_y, length_z );

  //===========================================================================
  // Shock Box Inputs
  //===========================================================================
#elif defined(SHOCK_BOX_2D)

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "shock_box_2d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 100;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 1.0;
  
  // the CFL and final solution time
  constexpr real_t CFL = 0.5;
  constexpr real_t final_time = 0.2;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[0] < 0.0 && x[1] < 0.0 ) {
        d = 0.125;
        p = 0.1;
      }
      else {
        d = 1.0;
        p = 1.0;
      }    
      return std::make_tuple( d, v, p );
    };
  
  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  //===========================================================================
  // Shock Box Inputs
  //===========================================================================
#elif defined(SHOCK_BOX_3D)

  using mesh_t = mesh_3d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "shock_box_3d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 100;

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 10;
  constexpr size_t num_cells_z = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 1.0;
  constexpr real_t length_z = 1.0;
  
  // the CFL and final solution time
  constexpr real_t CFL = 1/3.;
  constexpr real_t final_time = 0.2;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[0] < 0 && x[1] < 0 && x[2] < 0 ) {
        d = 0.125;
        p = 0.1;
      }
      else {
        d = 1.0;
        p = 1.0;
      }    
      return std::make_tuple( d, v, p );
    };
  
  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, num_cells_z, length_x, length_y, length_z );

#else

#  pragma message("NO CASE DEFINED!")

#endif

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  //===========================================================================
  // Parse arguments
  //===========================================================================

  // the usage stagement
  auto print_usage = [&argv]() {
    std::cout << "Usage: " << argv[0] 
              << " [--max_steps MAX_STEPS]"
              << " [--output_freq OUTPUT_FREQ]"
              << " [--case CASE_NAME]"
              << " [--postfix EXTENSION]"
              << std::endl << std::endl;
    std::cout << "\t--max_steps MAX_STEPS:\tOverride the number of "
              << "iterations with MAX_STEPS." << std::endl;
    std::cout << "\t--output_freq OUTPUT_FREQ:\tOverride the frequency "
              << "of output to occur ever OUTPUT_FREQ." << std::endl;
    std::cout << "\t--case CASE_NAME:\tOverride the case name "
              << "with CASE_NAME." << std::endl;
    std::cout << "\t--postfix EXTENSION:\tOverride the output extension "
              << "with EXTENSION." << std::endl;
  };

  while (1) {
    
    // Define the options
    static struct option long_options[] =
      {
        {"help",              no_argument, 0, 'h'},
        {"max_steps",   required_argument, 0, 'n'},
        {"output_freq", required_argument, 0, 'f'},
        {"case",        required_argument, 0, 'c'},
        {"postfix",     required_argument, 0, 'p'},
        {0, 0, 0, 0}
      };
      // getopt_long stores the option index here.
      int option_index = 0;

      auto c = getopt_long (argc, argv, "h", long_options, &option_index);

      // Detect the end of the options.
      if (c == -1) break;

      switch (c) {
        case 'c':
          short_warning ("Overriding case name with \"" << optarg << "\"" );
          prefix = optarg;
          break;

        case 'p':
          short_warning ("Overriding extension with \"" << optarg << "\"" );
          postfix = optarg;
          break;

        case 'f':
          short_warning ("Overriding output frequency to \"" << optarg << "\"" );
          output_freq = atoi( optarg );
          break;

        case 'n':
          short_warning ("Overriding max iterations to \"" << optarg << "\"" );
          max_steps = atoi( optarg );
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

  // success
  return 0;

}

