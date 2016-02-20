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
#include <iostream>
#include <utility>

// user includes
#include <ale/eqns/euler_eqns.h>
#include <ale/eos/ideal_gas.h>
#include <ale/geom/centroid.h>
#include <ale/geom/midpoint.h>
#include <ale/utils/zip.h>

#include <ale/mesh/burton/burton.h>
#include <ale/mesh/factory.h>

// some namespace aliases
namespace mesh  = ale::mesh;
namespace math  = ale::math;
namespace utils = ale::utils;
namespace geom  = ale::geom;


// math operations
using math::operator*;
using math::operator/;
using math::operator+;
using math::operator-;

// explicitly use some other stuff
using std::cout;
using std::cerr;
using std::endl;
using std::vector;


// mesh and some underlying data types
using mesh_t = mesh::burton_mesh_t;

using size_t = mesh_t::size_t;
using real_t = mesh_t::real_t;
using point_t = mesh_t::point_t;
using vector_t = mesh_t::vector_t;

using eos_t = ale::eos::ideal_gas_t<real_t>;

using eqns_t = ale::eqns::euler_eqns_t<real_t, mesh_t::dimension()>;
using state_data_t = eqns_t::state_data_t;
using flux_data_t = eqns_t::flux_data_t;

// hydro incdlues
#include "hydro_tasks.h"


///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

  //===========================================================================
  // Mesh Setup
  //===========================================================================

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 2;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.2;

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );
  cout << mesh;
  
  //===========================================================================
  // Field Creation
  //===========================================================================

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  register_state(mesh, "density",      cells,   real_t, persistent);
  register_state(mesh, "pressure",     cells,   real_t, persistent);
  register_state(mesh, "velocity",     cells, vector_t, persistent);

  register_state(mesh, "internal_energy", cells, real_t, persistent);
  register_state(mesh, "temperature",     cells, real_t, persistent);
  register_state(mesh, "sound_speed",     cells, real_t, persistent);

  // register state a global eos
  auto eos = register_global_state( mesh, "eos", eos_t );
  eos->set_gamma( 1.4 ); 
  eos->set_specific_heat_v( 1.0 ); 

  /*
  // register collections: the full independant+derived state
  register_collection( mesh, "full_state", 
    { "density", "velocity", "pressure", "internal_energy", "temperature", 
      "sound_speed" } );

  // register collections: the primitive (i.e. independant) state
  register_collection( mesh, "full_state", 
    { "density", "velocity", "pressure" } );
  */

  //===========================================================================
  // Initial conditions
  //===========================================================================

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
  
  // now call the main task to set the ics
  apps::hydro::initial_conditions( mesh, ics );
  
  //===========================================================================
  // Apply EOS
  //===========================================================================


  // now call the main task to set the ics
  apps::hydro::update_state_from_pressure( mesh );


  // now output the solution
  mesh::write_mesh("test/hydro_0.exo", mesh);

#if 0
  //===========================================================================
  // Flux Evaluation
  //===========================================================================

  // compute the fluxes.  here I am regestering a struct as the stored data
  // type since I will only ever be accesissing all the data at once.
  auto flux = register_state(mesh, "flux", edges, flux_data_t, temporary);

  // the flux function
  auto flux_function = [](const auto &wl, const auto &wr, const auto &n) { 
    auto fl = eqns_t::flux(wl, n);
    //auto fr = eqns_t::flux(wr, n);
    //auto f = fl + fr;
    return fl;
  };

  // TASK: loop over each edge and compute/store the flux
  // fluxes are stored on each edge

  for ( auto e : mesh.edges() ) {
    auto edge_id = e.id();

    // get the normal
    auto vs = mesh.vertices(e).to_vec();
    //auto normal = normal( vs[0]->coordinates(), vs[1]->coordinates() );
    vector_t normal;

    // get the cell neighbors
    auto c = mesh.cells(e).to_vec();
    auto cells = c.size();


    // get the left state
    auto left_cell = c[0];
    auto w_left = get_full_state( left_cell );    

    // interior cell
    if ( cells == 2 ) {
      auto right_cell = c[1];
      auto w_right = get_full_state( right_cell );
      flux[edge_id] = flux_function( w_left, w_right, normal );
    } 
    // boundary cell
    else {
      // make sure normal points outwards
      auto midpoint = geom::midpoint(e);
      auto centroid = geom::centroid(c[0]);
      auto delta = midpoint - centroid;
      //if ( dot_product( normal, delta ) < 0 ) normal = - normal;
      // compute the boundary flux
      //flux[edge_id] = eqns_t::flux( w_left, -normal );
    }
    

  }


  //===========================================================================
  // Residual evaluation
  //===========================================================================
    
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto cell : mesh.cells() ) {
    auto cell_id = cell.id();

    flux_data_t delta_u;
    math::fill(delta_u, 0.0);

    // loop over each connected edge
    for ( auto edge : mesh.edges(cell) ) {
      auto edge_id = edge.id();
      
      // get the cell neighbors
      auto neigh = mesh.cells(edge);
      auto num_neigh = neigh.size();
      auto neigh_it = neigh.begin();
      auto first_neigh_id= (*neigh_it).id();

      // figure out the sign
      int dir{1};
      if ( num_neigh == 2 ) 
        // interior
        dir = ( first_neigh_id == cell_id ) ? -1 : 1;
      else
        // boundary
        dir = 1;
      
      // add the contribution to this cell only
      math::plus_equal( delta_u, dir*flux[edge_id] );
      
    }

    

  }

  // now output the solution
  mesh::write_mesh("test/hydro_1.exo", mesh);

#endif

    
}



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
