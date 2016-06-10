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

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;

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
    auto vol = c->volume();
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

  // type aliases
  using eqns_t = eqns_t<T::num_dimensions>;

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

  // type aliases
  using eqns_t = eqns_t<T::num_dimensions>;
  using flux_data_t = flux_data_t<T::num_dimensions>;

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
template< typename T >
int32_t evaluate_time_step( T & mesh ) {

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;
  using eqns_t = eqns_t<T::num_dimensions>;
  using flux_data_t = flux_data_t<T::num_dimensions>;

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
    auto dVdt = eqns_t::volumetric_rate_of_change( dudt[c] );
    dti = std::abs(dVdt) / c->volume() / cfl->volume;
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
    cout << "accoustic";
    break;
  case 1:
    cout << "volume";
    break;
  case 2:
    cout << "growth";
    break;
  default:
    cout << "unknown";
    raise_runtime_error( "could not determine time step limit" );
    break;    
  };
  cout << endl;

  // invert dt and check against growth
  *time_step = *it;
      
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

  // type aliases
  using vector_t = typename T::vector_t;

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
template< typename T, typename BC >
int32_t evaluate_nodal_state( T & mesh, const BC & boundary_map ) {

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;
  using matrix_t = matrix_t< T::num_dimensions >; 
  using eqns_t = eqns_t<T::num_dimensions>;

  // get epsilon
  constexpr auto eps = std::numeric_limits<real_t>::epsilon();

  // get the number of dimensions and create a matrix
  constexpr size_t dims = T::num_dimensions;
  
  // access what we need
  auto cell_state = cell_state_accessor<T>( mesh );
  auto vertex_vel = access_state( mesh, "node_velocity", vector_t );
  auto Mpc = access_state( mesh, "corner_matrix", matrix_t );
  auto npc = access_state( mesh, "corner_normal", vector_t );


  // get the current time
  auto soln_time = mesh.time();

  //----------------------------------------------------------------------------
  // Loop over each vertex
  for ( auto vt : mesh.vertices() ) {

    // create the final matrix the point
    matrix_t Mp(0);
    vector_t rhs(0);

    //std::cout << vt.id() << std::endl;

    
    //--------------------------------------------------------------------------
    // build corner matrices
    for ( auto cn : mesh.corners(vt) ) {
      
      // initialize
      Mpc[cn] = 0;
      npc[cn] = 0;

      matrix_t Mc(0);
      vector_t nc(0);

      // a corner connects to a cell and a vertex
      auto cl = mesh.cells(cn).front();
      auto vt = mesh.vertices(cn).front();

      // get the cell state (there is only one)
      auto state = cell_state(cl);
      // the cell quantities
      auto dc = eqns_t::density( state );
      auto uc = eqns_t::velocity( state );
      auto pc = eqns_t::pressure( state );
      auto ac = eqns_t::sound_speed( state );
      auto Gc = eqns_t::impedance_multiplier( state );
    
      // get the nodal state
      auto uv = vertex_vel[vt];
      
      // iterate over the wedges in pairs
      auto wedges = mesh.wedges(cn);
      for ( auto wit = wedges.begin(); wit != wedges.end(); ++wit ) 
      {
        // get the first wedge normal
        auto n_plus  = (*wit)->facet_normal_right();
        auto l_plus  = abs(n_plus);
        assert( l_plus > 0 );
        auto un_plus  = n_plus / l_plus;
        
        bool symm_plus = false;
        matrix_t R_plus;
        if ( (*wit)->is_boundary() ) {
          auto f = mesh.faces(*wit).front();
          auto b = boundary_map.at( f->boundary_tag() );
          symm_plus = ( b->has_symmetry() );
          R_plus = math::reflection_matrix( un_plus );
        }
        
        ++wit;

        // get the second wedge normal
        auto n_minus = (*wit)->facet_normal_left();
        auto l_minus = abs(n_minus);
        assert( l_minus > 0 );
        auto un_minus = n_minus / l_minus;
        
        bool symm_minus = false;
        matrix_t R_minus;
        if ( (*wit)->is_boundary() ) {
          auto f = mesh.faces(*wit).front();
          auto b = boundary_map.at( f->boundary_tag() );
          symm_minus = ( b->has_symmetry() );
          R_minus = math::reflection_matrix( un_minus );
        }
        
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
        Mpc[cn] += M_plus;
        Mpc[cn] += M_minus;
        // compute the pressure coefficient
        npc[cn] += n_plus;
        npc[cn] += n_minus;
        // deal with symmetry ( but don't alter the true corner normal/matrices )

        ax_plus_y( M_plus, uc, rhs );      
        ax_plus_y( M_minus, uc, rhs );      

        if ( symm_plus ) {
          // nc -= n_plus;
          // math::matrix_multiply( M_plus,  R_plus,  Mc );
          Mp += M_plus;
          auto usymm = reflect( uc, n_plus );
          ax_plus_y( M_plus, usymm, rhs );      
        } else {
          rhs += pc * n_plus;
        }
        if ( symm_minus ) {
          // nc -= n_minus;
          // math::matrix_multiply( M_minus, R_minus, Mc );
          Mp += M_minus;
          auto usymm = reflect( uc, n_minus );
          ax_plus_y( M_minus, usymm, rhs );      
        } else {
          rhs += pc * n_minus;
        }

      } // wedges

      // add to the final corner matrix
      Mp += Mpc[cn];
    

      // add the pressure and velocity contributions to the system
      // rhs += pc * nc;
      // ax_plus_y( Mc, uc, rhs );      

    } // corners

    //--------------------------------------------------------------------------
    // now solve the system for the point velocity

    //---------- boundary point
    if ( vt->is_boundary() ) {

      // get the boundary tags
      const auto & point_tags =  vt->boundary_tags();

      // first check if this has a prescribed velocity.  If it does, then nothing to do
      auto vel_bc = std::find_if( 
        point_tags.begin(), point_tags.end(), 
        [&boundary_map](const auto & id) { 
          return boundary_map.at(id)->has_prescribed_velocity();
        } );
      if ( vel_bc != point_tags.end() ) {
        vertex_vel[vt] = boundary_map.at(*vel_bc)->velocity( vt->coordinates(), soln_time );
        continue;
      }

      // otherwise, apply the pressure conditions
      auto bnd_wedges = filter_boundary( mesh.wedges(vt) );
      for ( auto wit = bnd_wedges.begin(); wit != bnd_wedges.end(); ++wit ) 
      {
        // get the even wedge
        {
          auto f = mesh.faces(*wit).front();
          auto c = mesh.cells(*wit).front();
          auto b = boundary_map.at( f->boundary_tag() );
          if ( b->has_prescribed_pressure() ) {
            auto n = (*wit)->facet_normal_right();
            auto x = (*wit)->facet_centroid();
            rhs -= b->pressure( x, soln_time ) * n;
          }
        }
        // advance to the next wedge
        assert( wit != bnd_wedges.end() );
        ++wit;
        // get the odd wedge
        {
          auto f = mesh.faces(*wit).front();
          auto c = mesh.cells(*wit).front();
          auto b = boundary_map.at( f->boundary_tag() );
          if ( b->has_prescribed_pressure() ) {
            auto n = (*wit)->facet_normal_left();
            auto x = (*wit)->facet_centroid();
            rhs -= b->pressure( x, soln_time ) * n;
          }
        }
      }
      
    } // point type

    //---------- internal point
    // make sure sum(lpc) = 0
    //assert( abs(np) < eps && "error in norms" );
    // now solve for point velocity
    vertex_vel[vt] = math::solve( Mp, rhs );
  


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

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;
  using matrix_t = matrix_t< T::num_dimensions >; 
  using flux_data_t = flux_data_t<T::num_dimensions>;
  using eqns_t = eqns_t<T::num_dimensions>;

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
      auto pt = mesh.vertices(cn).front();
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

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;
  using flux_data_t = flux_data_t<T::num_dimensions>;
  using eqns_t = eqns_t<T::num_dimensions>;

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
    eqns_t::update_volume( u, cell->volume() );

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

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;

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

  // type aliases
  using vector_t = typename T::vector_t;

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

  // type aliases
  using vector_t = typename T::vector_t;

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

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;

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

  // type aliases
  using real_t = typename T::real_t;
  using vector_t = typename T::vector_t;

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

  auto cnt = mesh.time_step_counter();
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
