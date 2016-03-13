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
 * \file tasks.h
 * 
 * \brief Simple tasks related to solving full hydro solutions.
 *
 ******************************************************************************/

#pragma once

// hydro includes
#include "types.h"

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
  auto M = access_state( mesh, "cell:mass",       real_t );
  auto V = access_state( mesh, "cell:volume",     real_t );
  auto p = access_state( mesh, "cell:pressure",   real_t );
  auto v = access_state( mesh, "cell:velocity", vector_t );

  for ( auto c : mesh.cells() ) {
    // get the 
    auto x = c->centroid();
    // now copy the state to flexi
    real_t d;
    std::tie( d, v[c], p[c] ) = std::forward<F>(ics)( x );
    // set mass and volume now
    auto vol = c->area();
    M[c] = d*vol;
    V[c] = vol;
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t update_state_from_pressure( T & mesh ) {


  // get the collection accesor
  auto eos = access_global_state( mesh, "eos", eos_t );
  auto cell_state = cell_state_accessor<T>( mesh );

  for ( auto c : mesh.cells() ) {
    auto u = cell_state(c);
    eqns_t::update_state_from_pressure( u, *eos );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t update_state_from_energy( T & mesh ) {


  // get the collection accesor
  auto eos = access_global_state( mesh, "eos", eos_t );
  auto cell_state = cell_state_accessor<T>( mesh );

  for ( auto c : mesh.cells() ) {
    auto u = cell_state(c);
    eqns_t::update_state_from_energy( u, *eos );
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename T >
int32_t evaluate_time_step( T & mesh ) {

  // access what we need
  auto cell_state = cell_state_accessor<T>( mesh );

  auto dudt = access_state( mesh, "cell:residual", flux_data_t );

  auto delta_t = access_global_state( mesh, "time_step", real_t );
  real_t cfl = access_global_state( mesh, "cfl", real_t );
 
 
  //----------------------------------------------------------------------------
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_inv(0);

  for ( auto c : mesh.cells() ) {

    // get cell properties
    auto u = cell_state( c );
    auto lambda = E::fastest_wavespeed(u);

    // loop over each edge do find the min edge length
    auto edges = mesh.edges(c);
    assert( edges.size() > 2 && "not enough edges in cell" );
    auto eit = edges.begin();
    auto min_length = eit->length();
    std::for_each( ++eit, edges.end(), [&](auto && e) 
    { 
      min_length = std::min( e->length(), min_length );
    });

    // compute the inverse of the time scale
    auto dti =  lambda / min_length;
    // check for the maximum value
    dt_inv = std::max( dti, dt_inv );

    // now check the volume change
    auto dVdt = E::volumetric_rate_of_change( dudt[c] );
    dti = std::abs(dVdt) / c->area();
    // check for the maximum value
    dt_inv = std::max( dti, dt_inv );

  } // cell
  //----------------------------------------------------------------------------

  assert( dt_inv > 0 && "infinite delta t" );

  // invert dt and apply cfl
  delta_t = cfl / dt_inv;

  return 0;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute corner quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t evaluate_corner_coef( T & mesh ) {

  // get the number of dimensions and create a matrix
  constexpr size_t dims = T::dimension();
  
  // access what we need
  auto cell_state = cell_state_accessor<T>( mesh );

  auto Mpc = access_state( mesh, "corner:matrix", matrix_t );
  auto npc = access_state( mesh, "corner:normal", vector_t );

  //----------------------------------------------------------------------------
  // Loop over each corner
  for ( auto cn : mesh.corners() ) {

    // a corner connects to a cell and a vertex
    auto cl = cn->cell();

    // get the cell state (there is only one)
    auto state = cell_state(cl);
    // the cell quantities
    auto zc = eqns_t::impedance( state );
    
    // get the two wedge normals of the corner
    auto wedges = cn->wedges();
    auto n_plus  = wedges.front()->cell_facet_normal();
    auto n_minus = wedges.back() ->cell_facet_normal();

    // compute the mass matrix
    auto M_plus  = math::outer_product( n_plus , n_plus  );
    auto M_minus = math::outer_product( n_minus, n_minus );
    // the final matrix
    // Mpc = zc * ( lpc^- npc^-.npc^-  + lpc^+ npc^+.npc^+ );
    Mpc[cn] = M_plus + M_minus;
    Mpc[cn] *= zc;
    // compute the pressure coefficient
    npc[cn] = n_plus + n_minus;
  } // corner
  //----------------------------------------------------------------------------


  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t evaluate_nodal_state( T & mesh ) {

  // get epsilon
  constexpr auto eps = std::numeric_limits<real_t>::epsilon();

  // get the number of dimensions and create a matrix
  constexpr size_t dims = T::dimension();
  
  // access what we need
  auto cell_state = cell_state_accessor<T>( mesh );
  auto vertex_vel = access_state( mesh, "node:velocity", vector_t );
  auto Mpc = access_state( mesh, "corner:matrix", matrix_t );
  auto npc = access_state( mesh, "corner:normal", vector_t );

  //----------------------------------------------------------------------------
  // Loop over each vertex
  for ( auto vt : mesh.vertices() ) {

    // create the final matrix the point
    matrix_t Mp(0);
    vector_t rhs(0);
    vector_t np(0);

    //--------------------------------------------------------------------------
    // build point matrix
    for ( auto cn : mesh.corners(vt) ) {
      // corner attaches to one cell and one point
      auto cl = cn->cell();
      // get the cell state (there is only one)
      auto state = cell_state(cl);
      // the cell quantities
      auto pc = eqns_t::pressure( state );
      auto uc = eqns_t::velocity( state );
      // add to the global matrix
      Mp += Mpc[cn];
      // add to the normal vector
      np += npc[cn];
      // add the pressure and velocity contributions to the system
      rhs += pc * npc[cn];
      ax_plus_y( Mpc[cn], uc, rhs );      
    } // corner

    //--------------------------------------------------------------------------
    // now solve the system for the point velocity

    //---------- boundary point
    if ( vt->is_boundary() ) {
      
      // at boundary points, solve a larger system
      math::matrix<real_t, dims+1, dims+1> Mp_b;
      math::vector<real_t, dims+1> rhs_b;

      // create the new matrix
      // A = [ Mp       lpc*npc ]
      //     [ lpc*npc  0
      for ( auto i=0; i<dims; i++ ) {
        Mp_b( i, dims ) = np[i];
        Mp_b( dims, i ) = np[i];
        for ( auto j=0; j<dims; j++ )
          Mp_b( i, j ) = Mp(i,j);        
      }
      Mp_b( dims, dims ) = 0;

      // create the new rhs
      // y = [ pc*npc   bc ]
      std::copy( rhs.begin(), rhs.end(), rhs_b.begin() );
      rhs_b[ dims ] = 0;

      // now solve it
      auto res = math::solve( Mp_b, rhs_b );
      for ( auto i=0; i<dims; i++ ) vertex_vel[vt][i] = res[i];
      
    }
    //---------- internal point
    else {
      // make sure sum(lpc) = 0
      assert( abs(np) < eps && "error in norms" );
      // now solve for point velocity
      vertex_vel[vt] = math::solve( Mp, rhs );
    } // point type

  } // vertex
  //----------------------------------------------------------------------------


  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate residuals
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t evaluate_forces( T & mesh ) {

  // access what we need
  auto dudt = access_state( mesh, "cell:residual", flux_data_t );

  auto vertex_vel = access_state( mesh, "node:velocity", vector_t );

  auto Mpc = access_state( mesh, "corner:matrix", matrix_t );
  auto npc = access_state( mesh, "corner:normal", vector_t );

  auto cell_state = cell_state_accessor<T>( mesh );

  //----------------------------------------------------------------------------
  // TASK: loop over each cell and compute the residual
  for ( auto cl : mesh.cells() ) {

    // local cell residual
    math::fill( dudt[cl], 0.0 );

    // get the cell state (there is only one)
    auto state = cell_state(cl);
    // the cell quantities
    auto pc = eqns_t::pressure( state );
    auto uc = eqns_t::velocity( state );

    //--------------------------------------------------------------------------
    // compute subcell forces
    for ( auto cn : mesh.corners(cl) ) {

      // corner attaches to one point and zone
      auto pt = cn->vertex();
      // get the vertex velocity
      auto uv = vertex_vel[pt];

      // compute the subcell force
      // lpc*Pc*npc - Mpc*(uv-uc)
      auto force =  pc * npc[cn];
      auto delta_u = uc - uv;
      ax_plus_y( Mpc[cn], delta_u, force );      

      // the fluxes
      auto dmass =  math::dot_product( npc[cn], uv );
      auto dmom  = - force;
      auto dener = - math::dot_product( force, uv );

      // add contribution
      math::plus_equal( dudt[cl], flux_data_t{dmass, dmom, dener} );
      
    }// corners    
    //--------------------------------------------------------------------------
    
  } // cell
  //----------------------------------------------------------------------------

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t apply_update( T & mesh ) {

  // access what we need
  auto dudt = access_state( mesh, "cell:residual", flux_data_t );

  auto cell_state = cell_state_accessor<T>( mesh );

  // read only access
  real_t delta_t = access_global_state( mesh, "time_step", real_t );
 
  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto cell : mesh.cells() ) {

    // now compute the final update
    auto delta_u = delta_t * dudt[cell];

    // apply the update
    auto u = cell_state( cell );
    eqns_t::update_state_from_flux( u, delta_u );

  } // for
  //----------------------------------------------------------------------------

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to move the mesh
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t move_mesh( T & mesh ) {

  // access what we need
  auto vel = access_state( mesh, "node:velocity", vector_t );

  // read only access
  real_t delta_t = access_global_state( mesh, "time_step", real_t );
 
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto vt : mesh.vertices() )
    vt->coordinates() += delta_t * vel[vt];

  return 0;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief Output the solution
//!
//! \param [in] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t output( T & mesh, const std::string & prefix ) 
{
  static size_t cnt = 0;

  std::stringstream ss;
  ss << prefix;
  ss << std::setw( 7 ) << std::setfill( '0' ) << cnt++;
  ss << ".exo";
  
  mesh::write_mesh( ss.str(), mesh );
  
  return 0;
}

} // namespace hydro
} // namespace apps
