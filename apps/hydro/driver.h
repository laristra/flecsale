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
#include "tasks.h"
#include "types.h"
#include "../common/arguments.h"

// user includes
#include <flecsi/execution/reduction.h>
#include <ristra/utils/time_utils.h>
#include <ristra/io/catalyst/adaptor.h>


// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include<Kokkos_Core.hpp>

namespace apps {
namespace hydro {
  
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
  Kokkos::initialize();

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // get the client handle
  auto mesh = flecsi_get_client_handle(mesh_t, meshes, mesh0);

  // make sure geometry is up to date
  auto f = flecsi_execute_task(
	  update_geometry,
	  apps::hydro,
	  index,
	  mesh
	  );
  f.wait(); // DONT GO FORWARD UNTIL DONE!

  // get the input file
  auto args = apps::common::process_arguments( argc, argv );
  auto input_file_name =
    args.count("f") ? args.at("f") : std::string();

  // override any inputs if need be
  if ( !input_file_name.empty() ) {
    std::cout << "Using input file \"" << input_file_name << "\"."
              << std::endl;
    inputs_t::load( input_file_name );
  }

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
  if (!input_file_name.empty()) {
	  auto filename_char = flecsi_sp::utils::to_char_array(input_file_name);
	  f = flecsi_execute_task(
		  initial_conditions_from_file,
		  apps::hydro,
		  index,
		  mesh,
		  inputs_t::eos,
		  soln_time,
		  filename_char,
		  d, v, e, p, T, a);
  } else {
	  f = flecsi_execute_task(
		  initial_conditions,
		  apps::hydro,
		  index,
		  mesh,
		  inputs_t::eos,
		  soln_time,
		  d, v, e, p, T, a);
  }

  #ifdef HAVE_CATALYST
    auto insitu = io::catalyst::adaptor_t(catalyst_scripts);
    std::cout << "Catalyst on!" << std::endl;
  #endif

  //===========================================================================
  // Pre-processing
  //===========================================================================
  
  auto prefix_char = flecsi_sp::utils::to_char_array( inputs_t::prefix );
 	auto postfix_char =  flecsi_sp::utils::to_char_array( "exo" );

  // now output the solution
#if 0
  auto has_output = (inputs_t::output_freq > 0);
  if (has_output) {
    auto f = flecsi_execute_task(
      output,
 			apps::hydro,
 			index,
 			mesh,
 			prefix_char,
 			postfix_char,
			time_cnt,
      soln_time,
 			d, v, e, p, T, a
    );
    f.wait();
  }
#endif
  auto runtime = Legion::Runtime::get_runtime();
  auto ctx = Legion::Runtime::get_context();


  // dump connectivity
#if 0
  auto name = flecsi_sp::utils::to_char_array( inputs_t::prefix+".txt" );
  auto f = flecsi_execute_task(print, apps::hydro, index, mesh, name);
#endif

  f.wait();
  // start a clock
  auto tstart = ristra::utils::get_wall_time();

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  for ( 
    size_t num_steps = 0;
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   
runtime->begin_trace(ctx, 42);
    //-------------------------------------------------------------------------
    // compute the time step

    // we dont need the time step yet
    auto global_future_time_step = flecsi_execute_reduction_task(
      evaluate_time_step, apps::hydro, index, min, double, mesh, d, v, e, p, T,
			a, inputs_t::CFL, inputs_t::final_time - soln_time
    );

    //-------------------------------------------------------------------------
    // try a timestep

    // compute the fluxes
    flecsi_execute_task( evaluate_fluxes, apps::hydro, index, mesh,
        d, v, e, p, T, a, F );
 
    //auto time_step = global_future_time_step.get();

    // Loop over each cell, scattering the fluxes to the cell
    f = flecsi_execute_task( 
      apply_update, apps::hydro, index, mesh, inputs_t::eos,
      global_future_time_step, F, d, v, e, p, T, a
    );

runtime->end_trace(ctx, 42);
    //-------------------------------------------------------------------------
    // Post-process

    // update time
    //soln_time += time_step;
    time_cnt++;

    // output the time step
#if 0
    if ( rank == 0 ) {
      cout << std::string(80, '=') << endl;
      auto ss = cout.precision();
      cout.setf( std::ios::scientific );
      cout.precision(6);
      cout << "|  " << "Step:" << std::setw(10) << time_cnt
           << "  |  Time:" << std::setw(17) << soln_time 
     //      << "  |  Step Size:" << std::setw(17) << time_step
           << "  |" << std::endl;
      cout.unsetf( std::ios::scientific );
      cout.precision(ss);
    }
#endif

#ifdef HAVE_CATALYST
    if (!catalyst_scripts.empty()) {
      auto vtk_grid = mesh::to_vtk( mesh );
      insitu.process( 
        vtk_grid, soln_time, num_steps, (num_steps==inputs_t::max_steps-1)
      );
    }
#endif

#if 0
    // now output the solution
    if ( has_output && 
        (time_cnt % inputs_t::output_freq == 0 || 
         num_steps==inputs_t::max_steps-1 ||
         std::abs(soln_time-inputs_t::final_time) < epsilon
        )  
      ) 
    {
      flecsi_execute_task(
        output,
	 			apps::hydro,
 				index,
 				mesh,
	 			prefix_char,
 				postfix_char,
 				time_cnt,
        soln_time,
 				d, v, e, p, T, a
      );
    }

#endif
  }

  //===========================================================================
  // Post-process
  //===========================================================================
  f.wait();    
  auto tdelta = ristra::utils::get_wall_time() - tstart;

  if ( rank == 0 ) {

    cout << "Final solution time is " 
         << std::scientific << std::setprecision(2) << soln_time
         << " after " << time_cnt << " steps." << std::endl;

    
    std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
              << tdelta << "s." << std::endl;

  }

  // dump solution for verification
#if 0
  {
    //    auto name = flecsi_sp::utils::to_char_array( inputs_t::prefix+"-solution.txt" );
    //    flecsi_execute_task( dump, apps::hydro, index, mesh, time_cnt, soln_time,
    //        d, v, e, p, name );
  }
#endif

  // success if you reached here
  Kokkos::finalize();
  return 0;

}

} // namespace
} // namespace
