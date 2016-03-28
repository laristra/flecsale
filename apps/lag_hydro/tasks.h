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
  auto M = access_state( mesh, "cell_mass",       real_t );
  auto V = access_state( mesh, "cell_volume",     real_t );
  auto p = access_state( mesh, "cell_pressure",   real_t );
  auto v = access_state( mesh, "cell_velocity", vector_t );

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

  auto dudt = access_state( mesh, "cell_residual", flux_data_t );

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

  auto vertex_vel = access_state( mesh, "node_velocity", vector_t );

  auto dudt = access_state( mesh, "cell_residual", flux_data_t );

  auto time_step = access_global_state( mesh, "time_step", real_t );
  auto cfl = access_global_state( mesh, "cfl", time_constants_t );
 
 
  //----------------------------------------------------------------------------
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_acc_inv(0);
  real_t dt_vol_inv(0);

  for ( auto c : mesh.cells() ) {

    // get cell properties
    auto state = cell_state( c );
    auto uc = eqns_t::velocity( state );
    auto ac = eqns_t::sound_speed( state );
    auto Gc = eqns_t::impedance_multiplier( state );

    // initial wavespeed without extra term
    auto lambda = ac;
    // check for max wavespeed with extra impedance term
    for ( auto v : mesh.vertices(c) )
      lambda = std::max( ac + Gc*abs( vertex_vel[v] - uc ), lambda );

    // loop over each edge do find the min edge length
    auto min_length = c->min_length();

    // compute the inverse of the time scale
    auto dti =  lambda / min_length / cfl->accoustic;
    // check for the maximum value
    dt_acc_inv = std::max( dti, dt_acc_inv );

    // now check the volume change
    auto dVdt = E::volumetric_rate_of_change( dudt[c] );
    dti = std::abs(dVdt) / c->area() / cfl->volume;
    // check for the maximum value
    dt_vol_inv = std::max( dti, dt_vol_inv );

  } // cell
  //----------------------------------------------------------------------------


  assert( dt_acc_inv > 0 && dt_vol_inv > 0 && "infinite delta t" );
  
  // get the individual cfls
  auto dt_acc = cfl->accoustic / dt_acc_inv;
  auto dt_vol = cfl->volume / dt_vol_inv;
  auto dt_growth = cfl->growth*(*time_step);

  // find the minimum one
  auto dts = std::array<real_t,3> { dt_acc, dt_vol, dt_growth };
  auto it = std::min_element( dts.begin(), dts.end() );

  cout << "Time step limit: ";
  switch ( std::distance( dts.begin(), it ) ) {
  case 0:
    std::cout << "accoustic";
    break;
  case 1:
    std::cout << "volume";
    break;
  case 2:
    std::cout << "growth";
    break;
  default:
    std::cout << "unknown";
    raise_runtime_error( "could not determine time step limit" );
    break;    
  };
  cout << endl;

  // invert dt and check against growth
  *time_step = *it;
      
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
  constexpr size_t dims = T::num_dimensions();
  
  // access what we need
  auto cell_state = cell_state_accessor<T>( mesh );

  auto vertex_vel = access_state( mesh, "node_velocity", vector_t );

  auto Mpc = access_state( mesh, "corner_matrix", matrix_t );
  auto npc = access_state( mesh, "corner_normal", vector_t );

  //----------------------------------------------------------------------------
  // Loop over each corner
  for ( auto cn : mesh.corners() ) {

    // a corner connects to a cell and a vertex
    auto cl = cn->cell();
    auto vt = cn->vertex();

    // get the cell state (there is only one)
    auto state = cell_state(cl);
    // the cell quantities
    auto dc = eqns_t::density( state );
    auto uc = eqns_t::velocity( state );
    auto ac = eqns_t::sound_speed( state );
    auto Gc = eqns_t::impedance_multiplier( state );
    
    // get the nodal state
    auto uv = vertex_vel[vt];

    // get the two wedge normals of the corner
    auto wedges = cn->wedges();
    auto n_plus  = wedges.front()->facet_normal_right();
    auto n_minus = wedges.back() ->facet_normal_left();
    auto l_plus  = abs(n_plus);
    auto l_minus = abs(n_minus);
    auto un_plus  = unit(n_plus);
    auto un_minus = unit(n_minus);
       
    // estimate impedances
    auto delta_u = uv - uc;
    auto zc_minus = dc * ( ac + Gc * std::abs(dot_product(delta_u, un_minus)) );
    auto zc_plus  = dc * ( ac + Gc * std::abs(dot_product(delta_u, un_plus )) );

    // compute the mass matrix
    auto M_plus  = math::outer_product( un_plus , un_plus  );
    auto M_minus = math::outer_product( un_minus, un_minus );
    M_plus  *= zc_plus  * l_plus;
    M_minus *= zc_minus * l_minus;
    // the final matrix
    // Mpc = zc * ( lpc^- npc^-.npc^-  + lpc^+ npc^+.npc^+ );
    Mpc[cn] = M_plus + M_minus;
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
int32_t estimate_nodal_state( T & mesh ) {

  // access what we need
  auto cell_vel = access_state( mesh, "cell_velocity", vector_t );
  auto vertex_vel = access_state( mesh, "node_velocity", vector_t );

  //----------------------------------------------------------------------------
  // Loop over each vertex
  for ( auto v : mesh.vertices() ) {
    vertex_vel[v] = 0.;
    auto cells = mesh.cells(v);
    for ( auto c : cells ) vertex_vel[v] += cell_vel[c];
    vertex_vel[v] /= cells.size();
  } // vertex
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
  constexpr size_t dims = T::num_dimensions();
  
  // access what we need
  auto cell_state = cell_state_accessor<T>( mesh );
  auto vertex_vel = access_state( mesh, "node_velocity", vector_t );
  auto Mpc = access_state( mesh, "corner_matrix", matrix_t );
  auto npc = access_state( mesh, "corner_normal", vector_t );

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


#if 0
  //----------------------------------------------------------------------------
  // enforce the boundary conditions
  for ( auto vt : mesh.vertices() ) {
    for ( auto ed : mesh.edges(vt) ) {

      if ( ed->is_boundary() ) {

        // get the unit normal
        auto n = unit( ed->normal() );
        
        // force the normal direction to zero to be sure
        auto uv = vertex_vel[vt];
        auto uv_n = dot_product( uv, n ) * n;
        auto uv_t = uv - uv_n;

        // take out the normal component
        vertex_vel[vt] = uv_t;
      }

    } // edge
  } // vertex
  //----------------------------------------------------------------------------
#endif

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
  auto dudt = access_state( mesh, "cell_residual", flux_data_t );

  auto vertex_vel = access_state( mesh, "node_velocity", vector_t );

  auto Mpc = access_state( mesh, "corner_matrix", matrix_t );
  auto npc = access_state( mesh, "corner_normal", vector_t );

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
int32_t apply_update( T & mesh, real_t coef ) {

  // access what we need
  auto dudt = access_state( mesh, "cell_residual", flux_data_t );

  auto cell_state = cell_state_accessor<T>( mesh );

  // read only access
  real_t delta_t = access_global_state( mesh, "time_step", real_t );

  // the time step factor
  auto fact = coef * delta_t;

  real_t mass(0);
  vector_t mom(0);
  real_t ener(0);
 
  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto cell : mesh.cells() ) {

    // get the cell state
    auto u = cell_state( cell );
    
    // now compute the final update
    auto delta_u = fact * dudt[cell];

    // apply the update
    eqns_t::update_state_from_flux( u, delta_u );
    eqns_t::update_volume( u, cell->area() );

#ifdef DEBUG
    // post update sums
    auto vel = eqns_t::velocity(u);
    auto et = eqns_t::total_energy(u);
    auto m  = eqns_t::mass(u);
    mass += m;
    mom  += m * vel;
    ener += m * et;
#endif
    
  } // for
  //----------------------------------------------------------------------------

#ifdef DEBUG
  auto ss = cout.precision();
  cout.setf( std::ios::scientific );
  cout.precision(8);
  cout << "Sums: mass = " << mass << ", momentum = " << mom << ", energy = " << ener << endl;
  cout.unsetf( std::ios::scientific );
  cout.precision(ss);
#endif

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to move the mesh
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t move_mesh( T & mesh, real_t coef ) {

  // access what we need
  auto vel = access_state( mesh, "node_velocity", vector_t );

  // read only access
  real_t delta_t = access_global_state( mesh, "time_step", real_t );

  // the time step factor
  auto fact = coef * delta_t;
 
  // Loop over each cell, scattering the fluxes to the cell
  for ( auto vt : mesh.vertices() )
    vt->coordinates() += fact * vel[vt];

  return 0;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t save_coordinates( T & mesh ) {

  // access what we need
  auto coord0 = access_state( mesh, "node_coordinates", vector_t );

  // Loop over vertices
  for ( auto vt : mesh.vertices() )
    coord0[vt] = vt->coordinates();

  return 0;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t restore_coordinates( T & mesh ) {

  // access what we need
  auto coord0 = access_state( mesh, "node_coordinates", vector_t );

  // Loop over vertices
  for ( auto vt : mesh.vertices() )
    vt->coordinates() = coord0[vt];

  return 0;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t save_solution( T & mesh ) {

  // access what we need
  auto vel  = access_state( mesh, "cell_velocity",   vector_t );
  auto vel0 = access_state( mesh, "cell_velocity_0", vector_t );

  auto ener  = access_state( mesh, "cell_internal_energy",   real_t );
  auto ener0 = access_state( mesh, "cell_internal_energy_0", real_t );

  // Loop over cells
  for ( auto c : mesh.cells() ) {
    vel0[c] = vel[c];
    ener0[c] = ener[c];
  }

  return 0;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t restore_solution( T & mesh ) {

  // access what we need
  auto vel  = access_state( mesh, "cell_velocity",   vector_t );
  auto vel0 = access_state( mesh, "cell_velocity_0", vector_t );

  auto ener  = access_state( mesh, "cell_internal_energy",   real_t );
  auto ener0 = access_state( mesh, "cell_internal_energy_0", real_t );

  // Loop over cells
  for ( auto c : mesh.cells() ) {
    vel[c] = vel0[c];
    ener[c] = ener0[c];
  }

  return 0;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief Output the solution
//!
//! \param [in] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t output( T & mesh, 
                const std::string & prefix, 
                const std::string & postfix, 
                size_t output_freq ) 
{

  auto cnt = mesh.get_time_step_counter();
  if ( cnt % output_freq != 0 ) return 0;

  std::stringstream ss;
  ss << prefix;
  ss << std::setw( 7 ) << std::setfill( '0' ) << cnt++;
  ss << "."+postfix;
  
  mesh::write_mesh( ss.str(), mesh );
  
  return 0;
}

} // namespace hydro
} // namespace apps
