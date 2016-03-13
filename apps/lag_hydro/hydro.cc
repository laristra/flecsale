/*~-------------------------------------------------------------------------~~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
/*!
 *
 * \file hydro.cc
 * 
 * \brief Simple tests related to solving full hydro solutions.
 *
 ******************************************************************************/

// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

// user includes
#include <ale/mesh/factory.h>

// hydro incdludes
#include "tasks.h"
#include "types.h"

// everything is in the hydro namespace
using namespace apps::hydro;

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

  //===========================================================================
  // Inputs
  //===========================================================================

  // the case prefix
  std::string prefix = "hydro";

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 2;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr real_t CFL = 0.5;
  constexpr real_t final_time = 0.2;

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

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );
  cout << mesh;
  
  //===========================================================================
  // Field Creation
  //===========================================================================

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  register_state(mesh, "cell:mass",     cells,   real_t, persistent);
  register_state(mesh, "cell:volume",   cells,   real_t, persistent);
  register_state(mesh, "cell:pressure", cells,   real_t, persistent);
  register_state(mesh, "cell:velocity", cells, vector_t, persistent);

  register_state(mesh, "node:velocity", vertices, vector_t, persistent);

  register_state(mesh, "cell:density",         cells,   real_t, persistent);
  register_state(mesh, "cell:internal_energy", cells, real_t, persistent);
  register_state(mesh, "cell:temperature",     cells, real_t, persistent);
  register_state(mesh, "cell:sound_speed",     cells, real_t, persistent);

  // register the time step and set a cfl
  register_global_state( mesh, "time_step", real_t );
  register_global_state( mesh, "cfl", real_t ) = CFL;
  

  // register state a global eos
  register_global_state( mesh, "eos", eos_t ) = eos;

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
  apps::hydro::output(mesh, prefix);

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  // compute the fluxes.  here I am regestering a struct as the stored data
  // type since I will only ever be accesissing all the data at once.
  register_state(mesh, "cell:residual", cells, flux_data_t, temporary);
  register_state(mesh, "corner:matrix", corners, matrix_t, temporary);
  register_state(mesh, "corner:normal", corners, vector_t, temporary);

  // get the current time
  auto soln_time = mesh.get_time();

  // a counter for time steps
  auto time_cnt = 0;

  do {   

    // evaluate corner matrices and normals
    apps::hydro::evaluate_corner_coef( mesh );
    
    // compute the nodal velocity
    apps::hydro::evaluate_nodal_state( mesh );

    // compute the fluxes
    apps::hydro::evaluate_forces( mesh );

    // compute the time step
    apps::hydro::evaluate_time_step<eqns_t>( mesh );
    
    // access the computed time step and make sure its not too large
    auto time_step = access_global_state( mesh, "time_step", real_t );   
    time_step = std::min( *time_step, final_time - soln_time );       

    cout << "step =  " << std::setw(4) << time_cnt++
         << std::setprecision(2)
         << ", time = " << std::scientific << soln_time
         << ", dt = " << std::scientific << *time_step
         << std::endl;

    
    // Loop over each cell, scattering the fluxes to the cell
    apps::hydro::apply_update( mesh );

    // move the mesh
    apps::hydro::move_mesh( mesh );

    // Update derived solution quantities
    apps::hydro::update_state_from_energy( mesh );

    // update time
    soln_time += *time_step;
    mesh.set_time( soln_time );

    // now output the solution
    apps::hydro::output(mesh, prefix);

  } while ( soln_time < final_time );


  //===========================================================================
  // Post-process
  //===========================================================================
    

  // success
  return 0;

}



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
