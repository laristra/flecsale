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
#include <utility>

// user includes
#include <ale/common/types.h>
#include <ale/eqns/euler_eqns.h>
#include <ale/eos/ideal_gas.h>
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
using flexi::normal;
using flexi::persistent;
using flexi::temporary;

using namespace ale;



// mesh and some underlying data types
using size_t = std::size_t;
using real_t = common::real_t;
using vector_t = math::vector<real_t,2>;

using mesh_t = flexi::burton_mesh_t;

using eos_t = eos::ideal_gas_t<real_t>;

using eqns_t = eqns::euler_eqns_t<real_t, mesh_t::dimension()>;
using state_data_t = eqns_t::state_data_t;
using flux_data_t = eqns_t::flux_data_t;

using math::operator*;
using math::operator/;
using math::operator+;
using math::operator-;


///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
TEST(hydro, simple) {

  //===========================================================================
  // Mesh Setup
  //===========================================================================


  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 2;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.2;

  constexpr real_t x0 = -length_x/2.0;
  constexpr real_t y0 = -length_y/2.0;


  // this is the mesh object
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

  // register state also returns an accessor
  auto eos_data = register_state(mesh, "eos", cells, eos_t);


  // a lambda function to get the full independant+derived state
  auto get_full_state = [&mesh](auto cell_id) { 

    auto density  = access_state(mesh, "density",   real_t);
    auto pressure = access_state(mesh, "pressure",   real_t);
    auto velocity = access_state(mesh, "velocity", vector_t);
    auto energy      = access_state(mesh, "internal_energy", real_t);
    auto temperature = access_state(mesh, "temperature", real_t);
    auto sound_speed = access_state(mesh, "sound_speed", real_t);

    return std::forward_as_tuple( density [ cell_id ],
                                  velocity[ cell_id ],
                                  pressure[ cell_id ],
                                  energy  [ cell_id ],
                                  temperature[ cell_id ],
                                  sound_speed[ cell_id ] );
  };

  // a lambda function to get the primitive (i.e. independant) state
  auto get_prim_state = [&mesh](auto cell_id) { 

    auto density  = access_state(mesh, "density",   real_t);
    auto pressure = access_state(mesh, "pressure",   real_t);
    auto velocity = access_state(mesh, "velocity", vector_t);

    return std::forward_as_tuple( density [ cell_id ],
                                  velocity[ cell_id ],
                                  pressure[ cell_id ] );
  };



  //===========================================================================
  // Initial conditions
  //===========================================================================

  // set the intial conditions.  This is really a TASK.
  for ( auto c : mesh.cells() ) {
    auto cell_id = c.id();
    auto x = mesh.centroid(c);

    // velocity is constant for now
    vector_t v = {0.0, 0.0};

    // figure out the left and right density/pressure
    real_t d;
    real_t p;
    if ( x[0] < 0.0 ) {
      d = 1.0;
      p = 1.0;
    }
    else {
      d = 0.125;
      p = 0.1;
    }    

    // a constant eos is used for now
    eos_t eos( /* gamma */ 1.4, /* cv */ 1.0); 
    eos.set_ref_state_dp(d, 2.5);

    // now copy the state to flexi
    get_prim_state( cell_id ) = std::make_tuple( d, v, p );
    eos_data[cell_id] = eos;
  }

  //===========================================================================
  // Apply EOS
  //===========================================================================


  // TASK: apply the eos to the state
  for ( auto c : mesh.cells() ) {
    auto cell_id = c.id();

    auto u     = get_full_state(cell_id);
    auto & eos = eos_data [cell_id];

    eqns_t::update_state_from_pressure( u, eos );
  }



  // now output the solution
  flexi::write_mesh("test/hydro_0.exo", mesh);

  //===========================================================================
  // Flux Evaluation
  //===========================================================================

  // compute the fluxes.  here I am regestering a struct as the stored data
  // type since I will only ever be accesissing all the data at once.
  auto flux = register_state(mesh, "flux", edges, flux_data_t, temporary);

  // the flux function
  auto flux_function = [](const auto &wl, const auto &wr, const auto &n) { 
    auto fl = eqns_t::flux(wl, n);
    auto fr = eqns_t::flux(wr, n);
    auto f = 0.5 * ( fl + fr );
    return f;
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
    auto left_cell_id = c[0];
    auto w_left = get_full_state( left_cell_id );    

    // interior cell
    if ( cells == 2 ) {
      auto right_cell_id = c[1];
      auto w_right = get_full_state( right_cell_id );
      flux[edge_id] = flux_function( w_left, w_right, normal );
    } 
    // boundary cell
    else {
      // make sure normal points outwards
      auto midpoint = mesh.midpoint(e);
      auto centroid = mesh.centroid(c[0]);
      auto delta = midpoint - centroid;
      //if ( dot_product( normal, delta ) < 0 ) normal = - normal;
      // compute the boundary flux
      //flux[edge_id] = eqns_t::flux( w_left, -normal );
    }
    

  }

#if 0


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

      // figure out the sign
      int dir;
      auto it = neigh.begin();
      if ( num_neigh == 2 ) 
        dir = ( it->id() == cell_id ) ? -1 : 0;
      else
        dir = 1;
      
      // add the contribution to this cell only
      math::plus_equal( delta_u, dir*flux[edge_id] );
      
    }

    

  }

  // now output the solution
  flexi::write_mesh("test/hydro_1.exo", mesh);


#endif

    
}



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
