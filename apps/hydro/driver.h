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
  1, 
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro, 
  velocity,
  mesh_t::vector_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro,
  internal_energy,
  mesh_t::real_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t, 
  hydro, 
  pressure,
  mesh_t::real_t, 
  dense, 
  1, 
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,
  hydro,
  temperature,
  mesh_t::real_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,
  hydro,
  sound_speed,
  mesh_t::real_t,
  dense,
  1,
  mesh_t::index_spaces_t::cells
);

// Here I am regestering a struct as the stored data
// type since I will only ever be accesissing all the data at once.
flecsi_register_field(
  mesh_t, 
  hydro, 
  flux, 
  flux_data_t, 
  dense, 
  1,
  mesh_t::index_spaces_t::faces
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

  //===========================================================================
  // Access what we need
  //===========================================================================
  
  auto d  = flecsi_get_handle(mesh, hydro,  density,   real_t, dense, 0);
  //auto d0 = flecsi_get_handle(mesh, hydro,  density,   real_t, dense, 1);
  auto v  = flecsi_get_handle(mesh, hydro, velocity, vector_t, dense, 0);
  //auto v0 = flecsi_get_handle(mesh, hydro, velocity, vector_t, dense, 1);
  auto e  = flecsi_get_handle(mesh, hydro, internal_energy, real_t, dense, 0);
  //auto e0 = flecsi_get_handle(mesh, hydro, internal_energy, real_t, dense, 1);

  auto p  = flecsi_get_handle(mesh, hydro,        pressure,   real_t, dense, 0);
  auto T  = flecsi_get_handle(mesh, hydro,     temperature, real_t, dense, 0);
  auto a  = flecsi_get_handle(mesh, hydro,     sound_speed, real_t, dense, 0);

  auto F = flecsi_get_handle(mesh, hydro, flux, flux_data_t, dense, 0);

  //===========================================================================
  // Initial conditions
  //===========================================================================
 
  // the solution time starts at zero
  real_t soln_time{0};  
  size_t time_cnt{0}; 


  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  flecsi_execute_task( 
    initial_conditions, 
    single, 
    mesh, 
    inputs_t::ics,
    inputs_t::eos,
    soln_time,
    d, v, e, p, T, a
  );

  #ifdef HAVE_CATALYST
    auto insitu = io::catalyst::adaptor_t(catalyst_scripts);
    std::cout << "Catalyst on!" << std::endl;
  #endif

  //===========================================================================
  // Pre-processing
  //===========================================================================

  // now output the solution
  auto has_output = (inputs_t::output_freq > 0);
  if (has_output) {
    auto name = prefix + "_" + apps::common::zero_padded( 0 ) + ".exo";
    auto name_char = utils::to_trivial_string( name );
    auto f =
      flecsi_execute_task(output, single, mesh, name_char, d, v, e, p, T, a);
  }


  // dump connectivity
  auto name = utils::to_trivial_string( prefix+".txt" );
  auto f = flecsi_execute_task(print, single, mesh, name);
  f.wait();

  // start a clock
  auto tstart = utils::get_wall_time();

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  auto & min_reduction = 
    flecsi::execution::context_t::instance().min_reduction();

  for ( 
    size_t num_steps = 0;
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    //-------------------------------------------------------------------------
    // compute the time step

    auto local_future = flecsi_execute_task( 
      evaluate_time_step, single, mesh, d, v, e, p, T, a,
      inputs_t::CFL, inputs_t::final_time - soln_time
    );

    // FleCSI does not yet support Future handling
    // auto global_future =
    //	flecsi::execution::context_t::instance().reduce_min(local_future);
    //flecsi::execution::flecsi_future__<mesh_t::real_t> *flecsi_future_time_step =
    //    &global_future;
    mesh_t::real_t hacked_time_step = 1.0e-3;

    //-------------------------------------------------------------------------
    // try a timestep

    // compute the fluxes
    f = 
      flecsi_execute_task( evaluate_fluxes, single, mesh, d, v, e, p, T, a, F );
 
    // Loop over each cell, scattering the fluxes to the cell
    f = flecsi_execute_task( 
      //apply_update, single, mesh, inputs_t::eos, flecsi_future_time_step, F, d, v, e, p, T, a
      apply_update, single, mesh, inputs_t::eos, hacked_time_step, F, d, v, e, p, T, a
    );
    // Be triply-sure that timing doesn't start until everyone is timestepping
    if (num_steps == 0) {
      f.wait();
      tstart = utils::get_wall_time();
    }


    //-------------------------------------------------------------------------
    // Post-process

    // update time
    bool silence_warnings = true;  // All sub-tasks have been launched.  So, it is now
                                   // acceptable to wait on global reduction and we can
                                   // silence the warnings this would produce.
    size_t index = 0;
    // This is where future_step_step should be evaluated
    // mesh_t::real_t future_time_step = flecsi_future_time_step->get(index, silence_warnings);
    mesh_t::real_t future_time_step = hacked_time_step;

    soln_time += future_time_step;
    time_cnt++;

    // output the time step
    if ( rank == 0 ) {
      cout << std::string(80, '=') << endl;
      auto ss = cout.precision();
      cout.setf( std::ios::scientific );
      cout.precision(6);
      cout << "|  " << "Step:" << std::setw(10) << time_cnt+1
           << "  |  Time:" << std::setw(17) << soln_time + future_time_step
           << "  |  Step Size:" << std::setw(17) << future_time_step
           << "  |" << std::endl;
      cout.unsetf( std::ios::scientific );
      cout.precision(ss);
    }

#ifdef HAVE_CATALYST
    if (!catalyst_scripts.empty()) {
      auto vtk_grid = mesh::to_vtk( mesh );
      insitu.process( 
        vtk_grid, soln_time, num_steps, (num_steps==inputs_t::max_steps-1)
      );
    }
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
      auto f = flecsi_execute_task(output, single, mesh, name_char, d, v, e, p, T, a);
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
