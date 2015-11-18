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

  using size_t = common::size_t;
  using real_t = common::real_t;
  using vector_t = utils::vector_t<real_t,2>;

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 2;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.2;

  constexpr real_t x0 = -length_x/2.0;
  constexpr real_t y0 = -length_y/2.0;

  auto solution = [=](const auto &x, auto t) { 
    if ( x[0] < 0.0 ) 
      return std::make_tuple( real_t(1.0),   real_t(1.0), vector_t{0.0,0.0} ); 
    else 
      return std::make_tuple( real_t(0.125), real_t(0.1), vector_t{0.0,0.0} ); 
  };

  // mesh and some underlying data types
  using mesh_t = flexi::burton_mesh_t;

  using eqns_t = euler_eqns_t<mesh_t::dimensions>;
  using eos_t = eos::ideal_gas_t;


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


  // setup eos
  eos_t eos;
  auto get_pressure         = std::bind( &eos_t::compute_pressure,        std::cref(eos), _1, _2 );
  auto get_internal_energy  = std::bind( &eos_t::compute_internal_energy, std::cref(eos), _1, _2 );

  // set the intial conditions
  for ( auto c : mesh.cells() ) {
    auto cell_id = c->id();
    auto cell_center = mesh.centroid(c);

    std::tie( density [cell_id], 
              pressure[cell_id], 
              velocity[cell_id] ) = solution( cell_center, 0.0 );    
  }



  // compute the conserved quantities
  auto momentum     = register_state(mesh, "momentum",     cells, vector_t);
  auto total_energy = register_state(mesh, "total_energy", cells,   real_t);

  for ( auto c : mesh.cells() ) {
    auto cell_id = c->id();

    auto &d = density[cell_id];
    auto &p = pressure[cell_id];
    auto &v = velocity[cell_id];

    auto e = get_internal_energy( d, p );
    total_energy[cell_id] = e + 0.5 * dot( v, v );
    momentum[cell_id] = d * v;
  }


  // compute the fluxes
  auto     mass_flux = register_state(mesh,     "mass_flux", edges,   real_t, temporary);
  auto momentum_flux = register_state(mesh, "momentum_flux", edges, vector_t, temporary);
  auto   energy_flux = register_state(mesh,   "energy_flux", edges,   real_t, temporary);

  for ( auto e : mesh.edges() ) {
    auto edge_id = e->id();

    // get the cell neighbors
    auto c = mesh.cells(e).toVec();
    
    auto cell_id = c[0]->id();

    auto dl = density[cell_id];
    auto pl = pressure[cell_id];
    auto vl = velocity[cell_id];

    eqns_t::primitive_state_t left(dl,vl,pl);
    

  }



  // now output the solution
  std::string name("test/hydro.exo");
  flexi::write_mesh(name, mesh);

    
} // TEST_F



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
