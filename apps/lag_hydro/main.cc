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
//#define SODX
//#define SODY
//#define SODXY
#define SEDOV

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

  //===========================================================================
  // SOD-X Inputs
  //===========================================================================

#ifdef SODX

  // the case prefix
  std::string prefix = "sodx";
  std::string postfix = "exo";

  // output frequency
  constexpr size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;

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

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  //===========================================================================
  // SOD-Y Inputs
  //===========================================================================

#elif defined(SODY)

  // the case prefix
  std::string prefix = "sody";
  std::string postfix = "exo";

  // output frequency
  constexpr size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 100;

  constexpr real_t length_x = 0.1;
  constexpr real_t length_y = 1.0;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[1] < 0.0 ) {
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

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  //===========================================================================
  // SOD-XY Inputs
  //===========================================================================

#elif defined(SODXY)

  // the case prefix
  std::string prefix = "sodxy";
  std::string postfix = "exo";

  // output frequency
  constexpr size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;

  // rotate by 45 degrees
  constexpr real_t degrees = 45;

  // this is a lambda function to set the initial conditions
  auto ics = [&degrees] ( const auto & x )
    {
      // compute some angles
      auto radians = degrees * math::pi / 180;
      auto cos = std::cos( radians );
      auto sin = std::sin( radians );

      // compute rotation matrix
      vector_t xrot;
      xrot[0] = x[0]*cos - x[1]*sin;
      xrot[1] = x[0]*sin + x[1]*cos;

      // compute solution
      real_t d, p;
      vector_t v(0);
      if ( xrot[1] < 0.0 ) {
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

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );
  mesh::rotate( mesh, degrees );


  //===========================================================================
  // SEDOV Inputs
  //===========================================================================

#elif defined(SEDOV)

  // the case prefix
  std::string prefix  = "sedov";
  std::string postfix = "vtk";

  // output frequency
  constexpr size_t output_freq = 100;

  // the grid dimensions
  constexpr size_t num_cells_x = 30;
  constexpr size_t num_cells_y = 30;

  constexpr real_t length_x = 1.2;
  constexpr real_t length_y = 1.2;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 0.25, .volume = 0.1, .growth = 1.01 };
  constexpr real_t initial_time_step = 1.e-5;
  constexpr real_t final_time = 1.0;

  // the value of gamma
  constexpr real_t gamma = 1.4;

  // compute a reference volume
  auto vol = (length_x / num_cells_x) * (length_y / num_cells_y);

  // compute a radial size
  auto delta_r = 
    std::sqrt( math::sqr( length_x/num_cells_x ) + 
               math::sqr( length_y/num_cells_y ) );
  
  // this is a lambda function to set the initial conditions
  auto ics = [&delta_r,&vol] ( const auto & x )
    {
      constexpr real_t e0 = 0.244816;
      real_t d = 1.0;
      vector_t v = 0;
      real_t p = 1.e-6;
      auto r = sqrt( x[0]*x[0] + x[1]*x[1] );
      if ( r < delta_r  ) 
        p = (gamma - 1) * d * e0 / vol;
      return std::make_tuple( d, v, p );
    };
  
  // setup an equation of state
  eos_t eos(  gamma, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, 0.0, 0.0, length_x, length_y );

#else

#  pragma message("NO CASE DEFINED!")

#endif

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  cout << mesh;
  
  //===========================================================================
  // Field Creation
  //===========================================================================

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  register_state(mesh, "cell_mass",     cells,   real_t, persistent);
  register_state(mesh, "cell_volume",   cells,   real_t, persistent);
  register_state(mesh, "cell_pressure", cells,   real_t, persistent);
  register_state(mesh, "cell_velocity", cells, vector_t, persistent);

  register_state(mesh, "cell_density",         cells,   real_t, persistent);
  register_state(mesh, "cell_internal_energy", cells, real_t, persistent);
  register_state(mesh, "cell_temperature",     cells, real_t, persistent);
  register_state(mesh, "cell_sound_speed",     cells, real_t, persistent);

  register_state(mesh, "cell_residual", cells, flux_data_t, temporary);
  register_state(mesh, "cell_velocity_0", cells, vector_t, temporary);
  register_state(mesh, "cell_internal_energy_0", cells, real_t, temporary);

  register_state(mesh, "node_velocity", vertices, vector_t, persistent);
  register_state(mesh, "node_coordinates", vertices, vector_t, temporary);

  register_state(mesh, "corner_matrix", corners, matrix_t, temporary);
  register_state(mesh, "corner_normal", corners, vector_t, temporary);

  // register the time step and set a cfl
  register_global_state( mesh, "time_step", real_t ) = initial_time_step;
  register_global_state( mesh, "cfl", time_constants_t ) = CFL;
  

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

    //--------------------------------------------------------------------------
    // Begin Time step
    //--------------------------------------------------------------------------

    // Save solution at n=0
    apps::hydro::save_coordinates( mesh );
    apps::hydro::save_solution( mesh );

    //--------------------------------------------------------------------------
    // Predictor step : Evaluate Forces at n=0
    //--------------------------------------------------------------------------

    // estimate the nodal velocity at n=0
    apps::hydro::estimate_nodal_state( mesh );

    // evaluate corner matrices and normals at n=0
    apps::hydro::evaluate_corner_coef( mesh );
    
    // compute the nodal velocity at n=0
    apps::hydro::evaluate_nodal_state( mesh );

    // compute the fluxes
    apps::hydro::evaluate_forces( mesh );

    //--------------------------------------------------------------------------
    // Time step evaluation
    //--------------------------------------------------------------------------

    // compute the time step
    apps::hydro::evaluate_time_step<eqns_t>( mesh );
    
    // access the computed time step and make sure its not too large
    auto time_step = access_global_state( mesh, "time_step", real_t );   
    time_step = std::min( *time_step, final_time - soln_time );       

    auto ss = cout.precision();
    cout.setf( std::ios::scientific );
    cout.precision(2);
    cout << "step =  " << std::setw(4) << time_cnt+1
         << ", time = " << soln_time
         << ", dt = " << *time_step
         << std::endl;
    cout.unsetf( std::ios::scientific );
    cout.precision(ss);

#if 1 // set to 0 for first order

    //--------------------------------------------------------------------------
    // Move to n+1/2
    //--------------------------------------------------------------------------

    // move the mesh to n+1/2
    apps::hydro::move_mesh( mesh, 0.5 );
    
    // update solution to n+1/2
    apps::hydro::apply_update( mesh, 0.5 );

    // Update derived solution quantities
    apps::hydro::update_state_from_energy( mesh );

    //--------------------------------------------------------------------------
    // Corrector : Evaluate Forces at n=0
    //--------------------------------------------------------------------------

    // evaluate corner matrices and normals at n=0
    apps::hydro::evaluate_corner_coef( mesh );
    
    // compute the nodal velocity at n=0
    apps::hydro::evaluate_nodal_state( mesh );

    // compute the fluxes
    apps::hydro::evaluate_forces( mesh );

    //--------------------------------------------------------------------------
    // Move to n+1/2
    //--------------------------------------------------------------------------

    // restore the solution to n=0
    apps::hydro::restore_coordinates( mesh );
    apps::hydro::restore_solution( mesh );

#endif

    // move the mesh to n+1
    apps::hydro::move_mesh( mesh, 1.0 );
    
    // update solution to n+1
    apps::hydro::apply_update( mesh, 1.0 );

    // Update derived solution quantities
    apps::hydro::update_state_from_energy( mesh );

    //--------------------------------------------------------------------------
    // End Time step
    //--------------------------------------------------------------------------

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
