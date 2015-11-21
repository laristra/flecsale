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
#include <cinchtest.h>
#include <iostream>

// user includes
#include <ale/eqns/euler_eqns.h>
#include <ale/eos/ideal_gas.h>
#include <ale/utils/vector.h>
#include <ale/utils/zip.h>

#include <flexi/io/io.h>
#include <flexi/specializations/burton/burton.h>
#include <flexi/specializations/burton/burton_io_exodus.h>


// explicitly use some stuff
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

using flexi::dot;
using flexi::persistent;
using flexi::temporary;
    
using std::placeholders::_1;
using std::placeholders::_2;

using namespace ale;

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
TEST(hydro, simple) {

  //===========================================================================
  // Case Parameters
  //===========================================================================

  using size_t = std::size_t;
  using real_t = common::real_t;
  using vector_t = utils::vector_t<real_t,2>;

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 2;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.2;

  constexpr real_t x0 = -length_x/2.0;
  constexpr real_t y0 = -length_y/2.0;

  // mesh and some underlying data types
  using mesh_t = flexi::burton_mesh_t;

  using eqns_t = eqns::euler_eqns_t<real_t, mesh_t::dimension()>;
  using eos_t = eos::ideal_gas_t<real_t>;

  // set the initial conditions
  auto solution = [](const auto &x, auto t) { 
    eos_t eos( /* gamma */ 1.4, /* cv */ 1.0);
    real_t d, p;
    if ( x[0] < 0.0 ) {
      d = 1.0;
      p = 1.0;
    }
    else {
      d = 0.125;
      p = 0.1;
    }    
    eos.set_ref_state_dp(d, 2.5);
    return std::make_tuple( d, p, vector_t{0.0,0.0}, eos ); 
  };

  // the flux function
  auto flux_function = [](const auto &wl, const auto &wr, const auto &n) { 
    eqns_t::flux_data_t f;
    eqns_t::flux_data_t fl = wl.flux(n);
    eqns_t::flux_data_t fr = wr.flux(n);
    return fl;
  };



  //===========================================================================
  // Mesh Setup
  //===========================================================================

  mesh_t mesh;

  // reserve storage for the mesh
  auto num_vertex = ( num_cells_x + 1 ) * ( num_cells_y + 1 );
  mesh.init_parameters( num_vertex );
  
  
  // create the individual vertices
  using vertex_t = mesh_t::vertex_t;
  vector<vertex_t*> vs;
  
  auto delta_x = length_x / num_cells_x;
  auto delta_y = length_y / num_cells_y;

  auto num_vert_x = num_cells_x + 1;
  auto num_vert_y = num_cells_y + 1;

  for(size_t j = 0; j < num_vert_y; ++j) {
    auto y = y0 + j*delta_y;
    for(size_t i = 0; i < num_vert_x; ++i) {
      auto x = x0 + i*delta_x;
      auto v = mesh.create_vertex({x, y});
      vs.push_back(v);
    }
    
  }
  
  // define each cell
  auto index = [=](auto i, auto j) { return i + num_vert_x*j; };
  
  for(size_t j = 0; j < num_cells_y; ++j)
    for(size_t i = 0; i < num_cells_x; ++i) {
      auto c = 
        mesh.create_cell({
              vs[ index(i  , j  ) ],
              vs[ index(i+1, j  ) ],
              vs[ index(i+1, j+1) ],
              vs[ index(i  , j+1) ]});
    }
  
  
  // now finalize the mesh setup
  mesh.init();

  // create some field data
  auto density  = register_state(mesh, "density",      cells,   real_t, persistent);
  auto pressure = register_state(mesh, "pressure",     cells,   real_t, persistent);
  auto velocity = register_state(mesh, "velocity",     cells, vector_t, persistent);

  auto energy      = register_state(mesh, "internal_energy", cells, real_t, persistent);
  auto temperature = register_state(mesh, "temperature",     cells, real_t, persistent);
  auto sound_speed = register_state(mesh, "sound_speed",     cells, real_t, persistent);

  auto eos_data = register_state(mesh, "eos", cells, eos_t);

  //===========================================================================
  // Initial conditions
  //===========================================================================

  // set the intial conditions
  for ( auto c : mesh.cells() ) {
    auto cell_id = c->id();
    auto cell_center = mesh.centroid(c);


    auto & d = density [cell_id];
    auto & p = pressure[cell_id];
    auto & v = velocity[cell_id];
    auto & e = energy  [cell_id];
    auto & t = temperature[cell_id];
    auto & a = sound_speed[cell_id];
    auto & eos = eos_data [cell_id];

    std::tie( d, p, v, eos ) = solution( cell_center, 0.0 );    
    
    e = eos.compute_internal_energy_dp( d, p );
    t = eos.compute_temperature_de( d, e );
    a = eos.compute_sound_speed_de( d, e );

  }


  // now output the solution
  flexi::write_mesh("test/hydro_0.exo", mesh);


  //===========================================================================
  // Flux Evaluation
  //===========================================================================

  // a lambda function to get the state
  auto get_state = [&](auto cell_id) { 
    auto d = density [ cell_id ];
    auto p = pressure[ cell_id ];
    auto v = velocity[ cell_id ];
    auto e = energy  [ cell_id ];
    auto t = temperature[ cell_id ];
    auto a = sound_speed[ cell_id ];
    return eqns_t::state_t{d,v,p,e,t,a}; 
  };



  // compute the fluxes
  using flux_data_t = eqns_t::flux_data_t;
  auto flux = register_state(mesh, "flux", edges, flux_data_t, temporary);

  //---------------------------------------------------------------------------
  // loop over each edge and compute/store the flux
  // fluxes are stored on each edge

  for ( auto e : mesh.edges() ) {
    auto edge_id = e->id();

    // get the cell neighbors
    auto c = mesh.cells(e).toVec();
    auto cells = c.size();

    // get the left state
    auto left_cell_id = c[0]->id();
    auto w_left = get_state( left_cell_id );    

    // interior cell
    if ( cells == 2 ) {
      auto right_cell_id = c[1]->id();
      auto w_right = get_state( right_cell_id );
      vector_t normal;
      flux[edge_id] = flux_function( w_left, w_right, normal );
    } 
    // boundary cell
    else {
      vector_t normal;
      flux[edge_id] = w_left.flux( normal );
    }
    

  }

  //===========================================================================
  // Residual evaluation
  //===========================================================================
    
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto cell : mesh.cells() ) {
    auto cell_id = cell->id();

    flux_data_t delta_u(0.0);

    // loop over each connected edge
    for ( auto edge : mesh.edges(cell) ) {
      auto edge_id = edge->id();
      
      // get the cell neighbors
      auto neigh = mesh.cells(edge).toVec();
      auto num_neigh = neigh.size();

      // figure out the sign
      decltype(neigh) dir;
      if ( cells == 2 ) 
        dir = ( neigh[0]->id() == cell_id ) ? -1 : 0;
      else
        dir = 1;
      
      // add the contribution to this cell only
      delta_u += dir * flux[edge_id];
      
    }

    

  }

  // now output the solution
  flexi::write_mesh("test/hydro_1.exo", mesh);




    
}



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
