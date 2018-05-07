/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "globals.h"
#include "tasks.h"
#include "types.h"

// user includes
#include <ristra/utils/time_utils.h>


// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

namespace apps {
namespace hydro {
  
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

  std::cout << "STARTING DRIVER COMPUTATION!!!!" << std::endl;


  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // get the client handle 
  auto mesh = flecsi_get_client_handle(mesh_t, meshes, mesh0);
 
  // check the mesh
  flecsi_execute_task( 
    validate_mesh, 
    apps::hydro,
    single, 
    mesh
  );

  
  //===========================================================================
  // Some typedefs
  //===========================================================================

  using size_t = typename mesh_t::size_t;
  using real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t; 

  // get machine zero
  constexpr auto epsilon = std::numeric_limits<real_t>::epsilon();
  const auto machine_zero = std::sqrt(epsilon);

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
  
	auto uc0 = flecsi_get_handle(mesh, hydro, cell_velocity, vector_t, dense, 1);
  auto ec0 = flecsi_get_handle(mesh, hydro, cell_internal_energy, real_t, dense, 1);

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
  
  // install each boundary
  tag_t bc_key = 0;
  for ( const auto & bc_pair : inputs_t::bcs )
  {
    auto bc_type = bc_pair.first.get();
    auto bc_function = bc_pair.second; 
    flecsi_execute_task(
        install_boundary,
        apps::hydro,
        single,
        mesh,
        soln_time,
        bc_key++,
        bc_type,
        bc_function);
  }

  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  flecsi_execute_task( 
    initial_conditions, 
    apps::hydro,
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

  auto prefix_char = flecsi_sp::utils::to_char_array( inputs_t::prefix );
 	auto postfix_char =  flecsi_sp::utils::to_char_array( "exo" );

  // now output the solution
  auto has_output = (inputs_t::output_freq > 0);
  if (has_output) {
    flecsi_execute_task(
      output,
 			apps::hydro,
 			single,
 			mesh,
 			prefix_char,
 			postfix_char,
			time_cnt,
      soln_time,
 			dc, uc, ec, pc, Tc, ac
    );
  }



  // dump connectivity
  auto name = flecsi_sp::utils::to_char_array( inputs_t::prefix+".txt" );
  auto f = flecsi_execute_task(print, apps::hydro, single, mesh, name);
  f.wait();

  // start a clock
  auto tstart = ristra::utils::get_wall_time();

	// the initial time step
	auto time_step = inputs_t::initial_time_step;

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  for (
    size_t num_steps = 0;
    (num_steps < inputs_t::max_steps && soln_time < inputs_t::final_time); 
    ++num_steps 
  ) {   

    //--------------------------------------------------------------------------
    // Begin Time step
    //--------------------------------------------------------------------------

    // Save solution at n=0
    flecsi_execute_task( save_coordinates, apps::hydro, single, mesh, xn );
    flecsi_execute_task( save_solution, apps::hydro, single, mesh, uc, ec, uc0, ec0 );


    //--------------------------------------------------------------------------
    // Predictor step : Evaluate Forces at n=0
    //--------------------------------------------------------------------------

    // estimate the nodal velocity at n=0
    flecsi_execute_task(
			 estimate_nodal_state,
			 apps::hydro,
       single,
			 mesh, uc, un
		);

    // compute the nodal velocity at n=0
    flecsi_execute_task( 
      evaluate_nodal_state,
      apps::hydro,
      single,
      mesh,
      soln_time,
      Vc, Mc, uc, pc, dc, ec, Tc, ac,
      un, npc, Fpc
    );

    // compute the fluxes
    flecsi_execute_task(
			 evaluate_residual,
		 	 apps::hydro,
			 single,
			 mesh,
			 un, npc, Fpc, dUdt
     );

    //--------------------------------------------------------------------------
    // Time step evaluation
    //--------------------------------------------------------------------------

    // compute the time step
    auto local_time_step_future = flecsi_execute_task(
			evaluate_time_step,
	 		apps::hydro,
	 		single,
			mesh,
			inputs_t::CFL,
			time_step,
			ac, dUdt
 		);
    
    // now we need it
    time_step =
      flecsi::execution::context_t::instance().reduce_min(local_time_step_future);
    time_step = std::min( time_step, inputs_t::final_time - soln_time );       

		if ( rank == 0 ) {
      cout << std::string(44, '=') << endl;
      auto ss = cout.precision();
      cout.setf( std::ios::scientific );
      cout.precision(6);
      cout << "| " << std::setw(8) << "Step:"
           << " | " << std::setw(13) << "Time:"
           << " | " << std::setw(13) << "Step Size:"
           << " |" << std::endl;
      cout << "| " << std::setw(8) << time_cnt+1
           << " | " << std::setw(13) << soln_time + (time_step)
           << " | " << std::setw(13) << time_step
           << " |" << std::endl;
      cout.unsetf( std::ios::scientific );
      cout.precision(ss);
		}

// #define USE_FIRST_ORDER_TIME_STEPPING
#ifndef USE_FIRST_ORDER_TIME_STEPPING // set to 0 for first order

    //--------------------------------------------------------------------------
    // Move to n+1/2
		//--------------------------------------------------------------------------

    // move the mesh to n+1/2
    flecsi_execute_task(
			 move_mesh, 
 			 apps::hydro,
       single,
			 mesh, 
			 un,
			 0.5*time_step
     );

	 	// update solution to n+1/2
    flecsi_execute_task(
			 apply_update, 
 			 apps::hydro,
       single,
			 mesh, 
			 0.5*time_step,
			 dUdt,
       Vc, Mc, uc, pc, dc, ec, Tc, ac
     );

    // Update derived solution quantities
    flecsi_execute_task( 
      update_state_from_energy,
			apps::hydro,
			single,
			mesh,
			inputs_t::eos,
			Vc, Mc, uc, pc, dc, ec, Tc, ac 
		);

    //--------------------------------------------------------------------------
    // Corrector : Evaluate Forces at n=1/2
    //--------------------------------------------------------------------------

    // compute the nodal velocity at n=1/2
    flecsi_execute_task( 
      evaluate_nodal_state,
      apps::hydro,
      single,
      mesh,
      soln_time,
      Vc, Mc, uc, pc, dc, ec, Tc, ac,
      un, npc, Fpc
    );

    // compute the fluxes
    flecsi_execute_task(
			 evaluate_residual,
		 	 apps::hydro,
			 single,
			 mesh,
			 un, npc, Fpc, dUdt
     );

    //--------------------------------------------------------------------------
    // Move to n+1
    //--------------------------------------------------------------------------

	  // restore the solution to n=0
    flecsi_execute_task(
			restore_coordinates,
 			apps::hydro,
 			single,
 			mesh,
			xn
 		);

    flecsi_execute_task(
	 		restore_solution,
 			apps::hydro,
 			single,
 			mesh,
 			uc0, uc, ec0, ec
 		);

#endif // USE_FIRST_ORDER_TIME_STEPPING


    // move the mesh to n+1
    flecsi_execute_task(
			 move_mesh, 
 			 apps::hydro,
       single,
			 mesh, 
			 un,
			 time_step
     );
    
	 	// update solution to n+1/2
    flecsi_execute_task(
			 apply_update, 
 			 apps::hydro,
       single,
			 mesh, 
			 time_step,
			 dUdt,
       Vc, Mc, uc, pc, dc, ec, Tc, ac
     );

    // Update derived solution quantities
    flecsi_execute_task( 
      update_state_from_energy,
			apps::hydro,
			single,
			mesh,
			inputs_t::eos,
			Vc, Mc, uc, pc, dc, ec, Tc, ac 
		);


    //--------------------------------------------------------------------------
    // End Time step
    //--------------------------------------------------------------------------

    // update time
    soln_time += time_step;
    time_cnt++;
  
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
 				single,
 				mesh,
	 			prefix_char,
 				postfix_char,
 				time_cnt,
        soln_time,
 				dc, uc, ec, pc, Tc, ac
      );
    }

  } // for

  //===========================================================================
  // Post-process
  //===========================================================================

  auto tdelta = ristra::utils::get_wall_time() - tstart;

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
