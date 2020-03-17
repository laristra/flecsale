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

// Need this as we move towards libraryification
#include <flecsi-sp/burton/mesh_interface.h>

// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

namespace apps {
namespace hydro {
  
// create some field data.  Fields are registered as struct of arrays.
// this allows us to access the data in different patterns.
flecsi_register_color(hydro, solution_time, mesh_t::real_t, 1);

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
  // Mesh Setup
  //===========================================================================

  // get the client handle
  auto mesh = flecsi_get_client_handle(mesh_t, meshes, mesh0);

  // make sure geometry is up to date
  flecsi_execute_task(
	  update_geometry,
	  apps::hydro,
	  index,
	  mesh
	  ).wait(); // DONT GO FORWARD UNTIL DONE!

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
  real_t initial_soln_time{0};  
  size_t time_cnt{0}; 

  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  if (!input_file_name.empty()) {
	  auto filename_char = flecsi_sp::utils::to_char_array(input_file_name);
	  flecsi_execute_task(
		  initial_conditions_from_file,
		  apps::hydro,
		  index,
		  mesh,
		  inputs_t::eos,
		  initial_soln_time,
		  filename_char,
		  d, v, e, p, T, a);
  } else {
	  flecsi_execute_task(
		  initial_conditions,
		  apps::hydro,
		  index,
		  mesh,
		  inputs_t::eos,
		  initial_soln_time,
		  d, v, e, p, T, a);
  }

  #ifdef HAVE_CATALYST
    auto insitu = io::catalyst::adaptor_t(catalyst_scripts);
    std::cout << "Catalyst on!" << std::endl;
  #endif
    
  // keep track of current solution time with futures and a global handle
  // so that we do not block on the top-level-task

  auto solution_time_handle = flecsi_get_color(hydro, solution_time, mesh_t::real_t, 0);

  flecsi_execute_task( 
    init_soln_time, apps::hydro, index, initial_soln_time, solution_time_handle
  );

  //===========================================================================
  // Pre-processing
  //===========================================================================
  
  auto prefix_char = flecsi_sp::utils::to_char_array( inputs_t::prefix );
 	auto postfix_char =  flecsi_sp::utils::to_char_array( "exo" );

  // now output the solution
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
      			solution_time_handle,
 			d, v, e, p, T, a
    );
    f.wait();
  }


  // dump connectivity
  auto name = flecsi_sp::utils::to_char_array( inputs_t::prefix+".txt" );
  auto f = flecsi_execute_task(print, apps::hydro, index, mesh, name);
  f.wait();

  // start a clock
  auto tstart = ristra::utils::get_wall_time();

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  for ( 
    size_t num_steps = 0;
    (num_steps < inputs_t::max_steps
    // we cannot make this comparison without blocking on the top-level-task
    // && soln_time_handle < inputs_t::final_time
    ); 
    ++num_steps 
  ) {   
    //-------------------------------------------------------------------------
    // compute the time step

    // we dont need the time step yet
    auto global_future_time_step = flecsi_execute_reduction_task(
      evaluate_time_step, apps::hydro, index, min, double, mesh, d, v, e, p, T,
			a, inputs_t::CFL, inputs_t::final_time, solution_time_handle
    );

    //-------------------------------------------------------------------------
    // try a timestep

    // compute the fluxes
    flecsi_execute_task( evaluate_fluxes, apps::hydro, index, mesh,
        d, v, e, p, T, a, F );
 

    time_cnt++;

    // Loop over each cell, scattering the fluxes to the cell
    // update time
    f = flecsi_execute_task( 
      apply_update, apps::hydro, index, mesh, inputs_t::eos,
      global_future_time_step, F, d, v, e, p, T, a,
      time_cnt, initial_soln_time, solution_time_handle
    );

    //-------------------------------------------------------------------------
    // Post-process

#ifdef HAVE_CATALYST
    if (!catalyst_scripts.empty()) {
      auto vtk_grid = mesh::to_vtk( mesh );
      insitu.process( 
        vtk_grid, initial_soln_time, num_steps, (num_steps==inputs_t::max_steps-1)
      );
    }
#endif


    // now output the solution
    if ( has_output && 
        (time_cnt % inputs_t::output_freq == 0 || 
         num_steps==inputs_t::max_steps-1
	 // this condition cannot be tested unless we block in the top-level-task
	 // || std::abs(soln_time_handle-inputs_t::final_time) < epsilon
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
        			solution_time_handle,
 				d, v, e, p, T, a
      );
    }


  }

  //===========================================================================
  // Post-process
  //===========================================================================
  f.wait();    
  auto tdelta = ristra::utils::get_wall_time() - tstart;

  flecsi_execute_task( 
    print_soln_time, apps::hydro, index, tdelta, time_cnt, solution_time_handle);
  
  // output the solution if it wasn't on the last time step
  if ( has_output && 
      (time_cnt % inputs_t::output_freq != 0)  ) 
  {
    flecsi_execute_task(
      output,
      apps::hydro,
      index,
      mesh,
      prefix_char,
      postfix_char,
      time_cnt,
      solution_time_handle,
      d, v, e, p, T, a
    );
  }

  // dump solution for verification
  {
    auto name = flecsi_sp::utils::to_char_array( inputs_t::prefix+"-solution.txt" );
    flecsi_execute_task( dump, apps::hydro, index, mesh, time_cnt, solution_time_handle,
        d, v, e, p, name );
  }

  // success if you reached here
  return 0;

}

} // namespace
} // namespace
