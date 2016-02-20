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
 * \file hydro_tasks.h
 * 
 * \brief Simple tasks related to solving full hydro solutions.
 *
 ******************************************************************************/

#pragma once

namespace apps {
namespace hydro {

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename F >
int32_t initial_conditions( T & mesh, F && ics ) {

  // get the collection accesor
  //auto prim_state = access_collection( mesh, "prim_state" );

  for ( auto c : mesh.cells() ) {
    // get the 
    auto x = geom::centroid(c);
    // now copy the state to flexi
    //prim_state[c] = std::forward<F>(ics)( x );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t update_state_from_pressure( T & mesh ) {


  // get the collection accesor
  //auto full_state = access_collection( mesh, "full_state" );
  auto eos = access_global_state( mesh, "eos", eos_t );

  for ( auto c : mesh.cells() ) {
    //auto u = full_state[c];
    //eqns_t::update_state_from_pressure( u, eos );
  }

}

} // namespace hydro
} // namespace apps
