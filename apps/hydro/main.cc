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

// right now cases are hard coded
//#define SODX_2D
#define SODX_3D
//#define SHOCK_BOX_2D
//#define SHOCK_BOX_3D

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

  //===========================================================================
  // X-direction SOD Inputs
  //===========================================================================

#ifdef SODX_2D

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sodx";
  std::string postfix = "dat";

  // output frequency
  constexpr size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 1;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr real_t CFL = 1.0;
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
  std::string prefix = "sodx";
  std::string postfix = "dat";

  // output frequency
  constexpr size_t output_freq = 1;

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
  std::string prefix = "shock_box";
  std::string postfix = "plt";

  // output frequency
  constexpr size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 100;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 1.0;
  
  // the CFL and final solution time
  constexpr real_t CFL = 0.5;
  constexpr real_t final_time = 0.2;

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
  std::string prefix = "shock_box";
  std::string postfix = "exo";

  // output frequency
  constexpr size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 100;
  constexpr size_t num_cells_z = 100;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 1.0;
  constexpr real_t length_z = 1.0;
  
  // the CFL and final solution time
  constexpr real_t CFL = 0.5;
  constexpr real_t final_time = 0.2;

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
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, num_cells_z, length_x, length_y, length_z );

#else

#  pragma message("NO CASE DEFINED!")

#endif

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // this is the mesh object
  cout << mesh;
  
  //===========================================================================
  // Field Creation
  //===========================================================================

  // type aliases
  using eqns_t = eqns_t<mesh_t::num_dimensions()>;
  using flux_data_t = flux_data_t<mesh_t::num_dimensions()>;

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  register_state(mesh, "density",      cells,   real_t, persistent);
  register_state(mesh, "pressure",     cells,   real_t, persistent);
  register_state(mesh, "velocity",     cells, vector_t, persistent);

  register_state(mesh, "internal_energy", cells, real_t, persistent);
  register_state(mesh, "temperature",     cells, real_t, persistent);
  register_state(mesh, "sound_speed",     cells, real_t, persistent);


  // compute the fluxes.  here I am regestering a struct as the stored data
  // type since I will only ever be accesissing all the data at once.
  register_state(mesh, "flux", faces, flux_data_t, temporary);

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
  apps::hydro::output(mesh, prefix, postfix, 1);

  //===========================================================================
  // Residual Evaluation
  //===========================================================================

  // get the current time
  auto soln_time = mesh.time();
  auto time_cnt  = mesh.time_step_counter();

  do {   

    // compute the time step
    apps::hydro::evaluate_time_step<eqns_t>( mesh );
    
    // access the computed time step and make sure its not too large
    auto time_step = access_global_state( mesh, "time_step", real_t );   
    time_step = std::min( *time_step, final_time - soln_time );       

    cout << "step =  " << std::setw(4) << time_cnt+1
         << std::setprecision(2)
         << ", time = " << std::scientific << soln_time
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

  } while ( soln_time < final_time );


  //===========================================================================
  // Post-process
  //===========================================================================
    
  // now output the solution
  apps::hydro::output(mesh, prefix, postfix, 1);

  // success
  return 0;

}



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
