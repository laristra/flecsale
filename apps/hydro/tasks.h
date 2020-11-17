/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tasks related to solving full hydro solutions.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// hydro includes
#include "types.h"

// flecsi includes
#include <flecsi-sp/io/io_exodus.h>
#include <flecsi/execution/context.h>
#include <flecsi/execution/execution.h>
#include <ristra/utils/string_utils.h>

// system includes
#include <iomanip>

namespace apps {
namespace hydro {

////////////////////////////////////////////////////////////////////////////////
//! \brief Update mesh geometry
//!
//! \param [in] mesh the mesh object
////////////////////////////////////////////////////////////////////////////////
void update_geometry(
  client_handle_r<mesh_t> mesh
) {
	mesh.update_geometry();
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void initial_conditions(
  client_handle_r<mesh_t>  mesh,
  eos_t eos,
  real_t soln_time,
  dense_handle_w<real_t> d,
  dense_handle_w<vector_t> v,
  dense_handle_w<real_t> e,
  dense_handle_w<real_t> p,
  dense_handle_w<real_t> T,
  dense_handle_w<real_t> a
) {

  for ( auto c : mesh.cells( flecsi::owned ) ) {
    auto lid = c.id();
    std::tie( d(c), v(c), p(c) ) = inputs_t::initial_conditions(
      mesh, lid, soln_time );
    eqns_t::update_state_from_pressure(
      pack( c, d, v, p, e, T, a ),
      eos
    );
  }

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void initial_conditions_from_file(
  client_handle_r<mesh_t>  mesh,
  eos_t eos,
  real_t soln_time,
  char_array_t filename,
  dense_handle_w<real_t> d,
  dense_handle_w<vector_t> v,
  dense_handle_w<real_t> e,
  dense_handle_w<real_t> p,
  dense_handle_w<real_t> T,
  dense_handle_w<real_t> a
) {
	auto ics = inputs_t::get_initial_conditions(filename.str());

	// This doesn't work with lua input
	//#pragma omp parallel for
	for ( auto c : mesh.cells( flecsi::owned ) ) {
		std::tie( d(c), v(c), p(c) ) = ics( c->centroid(), soln_time );
		eqns_t::update_state_from_pressure(
			pack( c, d, v, p, e, T, a ),
			eos
			);
	}
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate gradients
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////

inline void gradient_2d(
    const matrix_t & dxdx,
    real_t denom,
    const vector_t & dudx,
    vector_t & grad )
{
  grad[0] = (dudx[0]*dxdx(1,1)-dudx[1]*dxdx(0,1)) * denom;
  grad[1] = (dudx[1]*dxdx(0,0)-dudx[0]*dxdx(0,1)) * denom;
}

inline void gradient_3d(
    const matrix_t & dxdx,
    real_t denom,
    real_t min1,
    real_t min2,
    real_t min3,
    const vector_t & dudx,
    vector_t & grad )
{
  grad[0] = (
         dudx[0]*min1
       - dxdx(0,1)*(dudx[1]*dxdx(2,2)-dxdx(1,2)*dudx[2])
       + dxdx(0,2)*(dudx[1]*dxdx(1,2)-dxdx(1,1)*dudx[2])
     ) * denom;
  grad[1] = (
         dxdx(0,0)*(dudx[1]*dxdx(2,2)-dxdx(1,2)*dudx[2])
       - dudx[0]*min2
       + dxdx(0,2)*(dudx[2]*dxdx(0,1)-dxdx(0,2)*dudx[1])
     ) * denom;
  grad[2] = (
         dxdx(0,0)*(dudx[2]*dxdx(1,1)-dxdx(1,2)*dudx[1])
       - dxdx(0,1)*(dudx[2]*dxdx(0,1)-dxdx(0,2)*dudx[1])
       + dudx[0]*min3
     ) * denom;
}

void evaluate_gradients( 
  client_handle_r<mesh_t> mesh,
  dense_handle_r<real_t> d,
  dense_handle_r<vector_t> v,
  dense_handle_r<real_t> e,
  dense_handle_r<real_t> p,
  dense_handle_r<real_t> a,
  dense_handle_w<vector_t> grad_d,
  dense_handle_w<array_of_vector_t> grad_v,
  dense_handle_w<vector_t> grad_e,
  dense_handle_w<vector_t> grad_p,
  dense_handle_w<vector_t> grad_a
) {
  
  constexpr int num_dims = mesh_t::num_dimensions;
  std::array<real_t, num_dims+4> u0, u1;
  std::array<real_t, num_dims+4> umin, umax;
  std::array<real_t, num_dims+4> phi;
  std::array<vector_t, num_dims+4> dudx, grad_u;
  matrix_t dxdx;
  constexpr real_t tol = 1.e-9;
  
  for ( auto c0 : mesh.cells( flecsi::owned ) ) {

    // the the current point and state
    const auto & x0 = c0->centroid();
    const auto & vel = v(c0);
    size_t i=0;
    u0[i++] = d(c0);
    for (int d=0; d<num_dims; ++d) u0[i++] = vel[d];
    u0[i++] = e(c0);
    u0[i++] = p(c0);
    u0[i++] = a(c0);

    for ( int i=0; i<u0.size(); ++i ) {
      umin[i] = u0[i];
      umax[i] = u0[i];
    }

    // initialize the temporary variables
    dxdx = 0;
    for (int i=0; i<dudx.size(); ++i)
      dudx[i] = 0;

    // get the neighbors
    auto neighbors = mesh.cell_vertex_neighbors(c0);

    // Main loop over neighbors
    for ( auto c1 : neighbors ) {

      // get the state and coordinates
      const auto & x1 = c1->centroid();
      const auto & vel = v(c1);
      size_t i=0;
      u1[i++] = d(c1);
      for (size_t d=0; d<num_dims; ++d) u1[i++] = vel[d];
      u1[i++] = e(c1);
      u1[i++] = p(c1);
      u1[i++] = a(c1);
    
      for ( int i=0; i<u1.size(); ++i ) {
        umin[i] = std::min( u1[i], umin[i] );
        umax[i] = std::max( u1[i], umax[i] );
      }

      // compute the terms necessary for gradients
      auto dx = x1 - x0;

      if constexpr (num_dims == 1) {
        dxdx(0,0) += dx[0]*dx[0];
      }
      else if constexpr (num_dims == 2) {
        dxdx(0,0) += dx[0]*dx[0];
        dxdx(0,1) += dx[0]*dx[1];
        dxdx(1,1) += dx[1]*dx[1];
      }
      else if constexpr (num_dims == 3) {
        dxdx(0,0) += dx[0]*dx[0];
        dxdx(0,1) += dx[0]*dx[1];
        dxdx(0,2) += dx[0]*dx[2];
        dxdx(1,1) += dx[1]*dx[1];
        dxdx(1,2) += dx[1]*dx[2];
        dxdx(2,2) += dx[2]*dx[2];
      }

      for ( int i=0; i<u0.size(); ++i ) {
        auto du = u1[i] - u0[i];
        for ( int j=0; j<num_dims; ++j )
          dudx[i][j] += du*dx[j];
      }

    } // neighbors

    // finish gradients
    auto & dddx = grad_d(c0);
    auto & dvdx = grad_v(c0);
    auto & dedx = grad_e(c0);
    auto & dpdx = grad_p(c0);
    auto & dadx = grad_a(c0);

    if constexpr (num_dims == 1) {
      auto denom = 1 / dxdx(0,0);
      for ( int j=0; j<u0.size(); ++j )
        grad_u[j][0] = dudx[j][0] * denom;
    }
    else if constexpr (num_dims == 2) {
      auto denom = 1 / ( dxdx(0,0)*dxdx(1,1)-dxdx(0,1)*dxdx(0,1) );
      for ( int j=0; j<u0.size(); ++j )
        gradient_2d( dxdx, denom, dudx[j], grad_u[j] );
    }
    else if constexpr (num_dims == 3) {
      // First compute minors
      auto min1 = dxdx(1,1)*dxdx(2,2)-dxdx(1,2)*dxdx(1,2);
      auto min2 = dxdx(0,1)*dxdx(2,2)-dxdx(1,2)*dxdx(0,2);
      auto min3 = dxdx(0,1)*dxdx(1,2)-dxdx(1,1)*dxdx(0,2);
      // Now determinants
      auto denom = 1 / (dxdx(0,0)*min1-dxdx(0,1)*min2+dxdx(0,2)*min3);
      // compute
      for ( int j=0; j<u0.size(); ++j )
        gradient_3d( dxdx, denom, min1, min2, min3, dudx[j], grad_u[j] );
    }
    
    // limiters
    for ( int j=0; j<u0.size(); ++j ) phi[j] = 1;

    for ( auto v : mesh.vertices(c0) ) {
      
      // initialize vertex state
      const auto & xv = v->coordinates();
      for ( int j=0; j<u0.size(); ++j ) u1[j] = u0[j];
     
      // interpolate
      for (int i=0; i<num_dims; ++i ) {
        auto dx = xv[i] - x0[i];
        for ( int j=0; j<u1.size(); ++j ) u1[j] += dx*grad_u[j][i];
      }

      // compute limit
      for ( int j=0; j<u1.size(); ++j ) {
        auto du = u1[j] - u0[j];
        if ( du > tol ) {
          phi[j] = std::min(phi[j], (umax[j]-u0[j])/du);
        }
        else if ( du < -tol ) {
          phi[j] = std::min(phi[j], (umin[j]-u0[j])/du);
        }
      }

    } // verts

    // final gradient
    for ( int j=0; j<u0.size(); ++j ) {
      phi[j] = std::max<real_t>( phi[j], 0 );
      grad_u[j] *= phi[j];
    }

    // copy
    int j = 0;
    
    dddx = grad_u[j];
    ++j;
    
    for (int d=0; d<num_dims; ++d) {
      dvdx[d] = grad_u[j];
      ++j;
    }
    
    dedx = grad_u[j];
    ++j;
    
    dpdx = grad_u[j];
    ++j;
    
    dadx = grad_u[j];
  }

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size.
//!
//! \tparam E  The equation of state object to use.
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
real_t evaluate_time_step(
  client_handle_r<mesh_t> mesh,
  dense_handle_r<real_t> d,
  dense_handle_r<vector_t> v,
  dense_handle_r<real_t> e,
  dense_handle_r<real_t> p,
  dense_handle_r<real_t> T,
  dense_handle_r<real_t> a,
  real_t CFL,
  real_t final_time,
  color_handle_r<real_t> solution_time
) {
  real_t max_dt = final_time - solution_time;
 
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_inv(0);

  for ( auto c : mesh.cells( flecsi::owned ) ) {

    // get the solution state
    auto u = pack( c, d, v, p, e, T, a );

    // loop over each face
    for ( auto f : mesh.faces(c) ) {
      // estimate the length scale normal to the face
      auto delta_x = c->volume() / f->area();
      // compute the inverse of the time scale
      auto dti = eqns_t::fastest_wavespeed( u, f->normal() ) / delta_x;
      // check for the maximum value
      dt_inv = std::max( dti, dt_inv );
    } // edge

  } // cell

  if ( dt_inv <= 0 ) 
    THROW_RUNTIME_ERROR( "infinite delta t" );

  real_t time_step = 1 / dt_inv;
  time_step *= CFL;

  // access the computed time step and make sure its not too large
  time_step = std::min( time_step, max_dt );

  return time_step;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate fluxes at each face.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void evaluate_fluxes( 
  client_handle_r<mesh_t> mesh,
  dense_handle_r<real_t> d,
  dense_handle_r<vector_t> v,
  dense_handle_r<real_t> e,
  dense_handle_r<real_t> p,
  dense_handle_r<real_t> T,
  dense_handle_r<real_t> a,
  dense_handle_r<vector_t> grad_d,
  dense_handle_r<array_of_vector_t> grad_v,
  dense_handle_r<vector_t> grad_e,
  dense_handle_r<vector_t> grad_p,
  dense_handle_r<vector_t> grad_a,
  dense_handle_w_all<flux_data_t> flux
) {

  using namespace ristra::math;

  constexpr int num_dims = mesh_t::num_dimensions;
  const auto & face_list = mesh.faces( flecsi::owned );
  auto num_faces = face_list.size();

  for ( counter_t fit = 0; fit < num_faces; ++fit )
  {

    const auto & f = face_list[fit];
    const auto & fx = f->centroid();
    
    // get the cell neighbors
    const auto & cells = mesh.cells(f);
    auto num_cells = cells.size();

    // get the left state
    const auto & c0 = cells[0];
    const auto & cx0 = c0->centroid();
    auto w_left = pack_copy( c0, d, v, p, e, T, a );

    // reconstruct
    auto dx = fx - cx0;

    for ( int i=0; i<num_dims; ++i ) {
      std::get<0>(w_left) += dx[i]*grad_d(c0)[i];
      for ( int d=0; d<num_dims; ++d ) 
        std::get<1>(w_left)[d] += dx[i]*grad_v(c0)[d][i];
      std::get<2>(w_left) += dx[i]*grad_e(c0)[i];
      std::get<3>(w_left) += dx[i]*grad_p(c0)[i];
      std::get<5>(w_left) += dx[i]*grad_a(c0)[i];
    }
    
    // compute the face flux
    //
    // interior cell
    if ( num_cells == 2 ) {
      // get right state
      const auto & c1 = cells[1];
      const auto & cx1 = c1->centroid();
      auto w_right = pack_copy( c1, d, v, p, e, T, a );
      // reconstruct
      dx = fx - cx1;
      for ( int i=0; i<num_dims; ++i ) {
        std::get<0>(w_right) += dx[i]*grad_d(c1)[i];
        for ( int d=0; d<num_dims; ++d ) 
          std::get<1>(w_right)[d] += dx[i]*grad_v(c1)[d][i];
        std::get<2>(w_right) += dx[i]*grad_e(c1)[i];
        std::get<3>(w_right) += dx[i]*grad_p(c1)[i];
        std::get<5>(w_right) += dx[i]*grad_a(c1)[i];
      }
      // evaluate flux
      flux(f) = flux_function<eqns_t>( w_left, w_right, f->normal() );
    } 
    // boundary cell
    else {
      flux(f) = boundary_flux<eqns_t>( w_left, f->normal() );
    }

    // scale the flux by the face area
    flux(f) *= f->area();


  } // for
  //----------------------------------------------------------------------------

  // calculate ghost faces needed by shared cells
  // while duplicating shared faces calculations
  const auto & cell_list = mesh.cells( flecsi::shared );
  auto ncells = cell_list.size();

  #pragma omp parallel for
  for ( counter_t cit = 0; cit < ncells; ++cit )
  {

    const auto & c = cell_list[cit];

    // loop over each connected edge
    for ( auto f : mesh.faces(c) ) {
    
    const auto & fx = f->centroid();
    
    // get the cell neighbors
    const auto & cells = mesh.cells(f);
    auto num_cells = cells.size();

    // get the left state
    const auto & c0 = cells[0];
    const auto & cx0 = c0->centroid();
    auto w_left = pack_copy( c0, d, v, p, e, T, a );

    // reconstruct
    auto dx = fx - cx0;

    for ( int i=0; i<num_dims; ++i ) {
      std::get<0>(w_left) += dx[i]*grad_d(c0)[i];
      for ( int d=0; d<num_dims; ++d ) 
        std::get<1>(w_left)[d] += dx[i]*grad_v(c0)[d][i];
      std::get<2>(w_left) += dx[i]*grad_e(c0)[i];
      std::get<3>(w_left) += dx[i]*grad_p(c0)[i];
      std::get<5>(w_left) += dx[i]*grad_a(c0)[i];
    }
    
    // compute the face flux
    //
    // interior cell
    if ( num_cells == 2 ) {
      // right state
      const auto & c1 = cells[1];
      const auto & cx1 = c1->centroid();
      auto w_right = pack_copy( c1, d, v, p, e, T, a );
      // reconstruct
      dx = fx - cx1;
      for ( int i=0; i<num_dims; ++i ) {
        std::get<0>(w_right) += dx[i]*grad_d(c1)[i];
        for ( int d=0; d<num_dims; ++d ) 
          std::get<1>(w_right)[d] += dx[i]*grad_v(c1)[d][i];
        std::get<2>(w_right) += dx[i]*grad_e(c1)[i];
        std::get<3>(w_right) += dx[i]*grad_p(c1)[i];
        std::get<5>(w_right) += dx[i]*grad_a(c1)[i];
      }
      flux(f) = flux_function<eqns_t>( w_left, w_right, f->normal() );
    } 
    // boundary cell
    else {
      flux(f) = boundary_flux<eqns_t>( w_left, f->normal() );
    }
   
    // scale the flux by the face area
    flux(f) *= f->area();

    } // face

  } // for shared cell
  //----------------------------------------------------------------------------
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Keep a running future of the solution time for final output
//!
//! \param [in] initial solution time
//! \param [in] global time step update
//! \return future of current solution time
////////////////////////////////////////////////////////////////////////////////
real_t update_soln_time( 
  size_t time_cnt,
  real_t initial_soln_time,
  real_t delta_t,
  color_handle_rw<real_t> color_soln_time
) {

  real_t old_time = color_soln_time;

  real_t soln_time = initial_soln_time + old_time + delta_t;

  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  if (rank == 0) {
    // output the time step
    cout << std::string(80, '=') << endl;
    auto ss = cout.precision();
    cout.setf( std::ios::scientific );
    cout.precision(6);
    cout << "|  " << "Step:" << std::setw(10) << time_cnt
        << "  |  Time:" << std::setw(17) << soln_time 
        << "  |  Step Size:" << std::setw(17) << delta_t
        << "  |" << std::endl;
    cout.unsetf( std::ios::scientific );
    cout.precision(ss);
  }
 
  color_soln_time = soln_time;

  return soln_time;
}


template<typename T>
using handle_t =
  flecsi::execution::flecsi_future<T, flecsi::execution::launch_type_t::single>;

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution in each cell.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void apply_update( 
  client_handle_r<mesh_t> mesh,
  eos_t eos,
  handle_t<real_t> future_delta_t,
  dense_handle_r<flux_data_t> flux,
  dense_handle_rw<real_t> d,
  dense_handle_rw<vector_t> v,
  dense_handle_rw<real_t> e,
  dense_handle_rw<real_t> p,
  flecsi::dense_accessor<real_t, flecsi::rw, flecsi::rw, flecsi::wo> T,  // ghost never used
  dense_handle_rw<real_t> a,
  size_t time_cnt,
  real_t initial_soln_time,
  color_handle_rw<real_t> color_soln_time,
  real_t factor
) {

  //----------------------------------------------------------------------------
  // Loop over each cell, scattering the fluxes to the cell


  //auto delta_t = static_cast<real_t>( time_step );
  real_t delta_t = future_delta_t;
  delta_t*= factor;

  update_soln_time(time_cnt, initial_soln_time, delta_t, color_soln_time);

  const auto & cell_list = mesh.cells( flecsi::owned );
  auto num_cells = cell_list.size();

  #pragma omp parallel for
  for ( counter_t cit = 0; cit < num_cells; ++cit )
  {

    const auto & c = cell_list[cit];

    // initialize the update
    flux_data_t delta_u( 0 );

    // loop over each connected edge
    for ( auto f : mesh.faces(c) ) {
      
      // get the cell neighbors
      auto neigh = mesh.cells(f);
      auto num_neigh = neigh.size();

      // add the contribution to this cell only
      if ( neigh[0] == c )
        delta_u -= flux(f);
      else
        delta_u += flux(f);

    } // edge

    // now compute the final update
    delta_u *= delta_t/c->volume();

    // apply the update
    auto u = pack(c, d, v, p, e, T, a);
    eqns_t::update_state_from_flux( u, delta_u );

    // update the rest of the quantities
    eqns_t::update_state_from_energy( u, eos );

    // check the solution quantities
    if ( eqns_t::internal_energy(u) < 0 || eqns_t::density(u) < 0 ) 
      THROW_RUNTIME_ERROR( "Negative density or internal energy encountered!" );

  } // for
  //----------------------------------------------------------------------------
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Initialize future of the solution time for final output
//!
//! \param [in] initial solution time
//! \return future of current solution time
////////////////////////////////////////////////////////////////////////////////
void init_soln_time( 
  real_t initial_soln_time,
  color_handle_rw<real_t> soln_time
) {

  soln_time = initial_soln_time;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief Print final solution time
//!
//! \param [in] wall clock time
//! \param [in] number of time steps
//! \param [in] solution time
////////////////////////////////////////////////////////////////////////////////
void print_soln_time( 
  real_t tdelta,
  size_t time_cnt,
  color_handle_r<real_t> soln_time
) {
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  if (rank == 0) {
    cout << "Final solution time is " 
       << std::scientific << std::setprecision(2) << soln_time
       << " after " << time_cnt << " steps." << std::endl;

    std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
            << tdelta << "s." << std::endl;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
////////////////////////////////////////////////////////////////////////////////
void output( 
  client_handle_r<mesh_t> mesh, 
  char_array_t prefix,
	char_array_t postfix,
	size_t iteration,
  color_handle_r<real_t> time,
  dense_handle_r<real_t> d,
  dense_handle_r<vector_t> v,
  dense_handle_r<real_t> e,
  dense_handle_r<real_t> p,
  dense_handle_r<real_t> T,
  dense_handle_r<real_t> a
) {
  clog(info) << "OUTPUT MESH TASK" << std::endl;
 
  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // figure out this ranks file name
  auto output_filename = 
    prefix.str() + "_rank" + apps::common::zero_padded(rank) +
    "." + apps::common::zero_padded(iteration) + "." + postfix.str();

  // now outut the mesh
  using field_type = decltype(d);
  std::vector<field_type*> var_ptrs{&d, &p};
  std::vector<std::string> var_names{"density", "pressure"};
  flecsi_sp::io::io_exodus<mesh_t>::write(
	  output_filename, mesh, time, var_ptrs, var_names
  );
}

////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
////////////////////////////////////////////////////////////////////////////////
void print( 
  client_handle_r<mesh_t> mesh,
  char_array_t filename
) {

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  clog(info) << "PRINT MESH ON RANK " << rank << std::endl;
 
  // figure out this ranks file name
  auto name_and_ext = ristra::utils::split_extension( filename.str() );
  auto output_filename = 
    name_and_ext.first + "_rank" + apps::common::zero_padded(rank) +
    "." + name_and_ext.second;

  // dump to file
  std::cout << "Dumping connectivity to: " << output_filename << std::endl;
  std::ofstream file( output_filename );
  mesh.dump( file );

  // close file
  file.close();
  
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Dump solution to file for regression testing
////////////////////////////////////////////////////////////////////////////////
void dump(
	client_handle_r<mesh_t> mesh,
	size_t iteration,
  color_handle_r<real_t> future_time,
	dense_handle_r<real_t> d,
	dense_handle_r<vector_t> v,
	dense_handle_r<real_t> e,
	dense_handle_r<real_t> p,
	char_array_t filename) {

  real_t time = future_time;

	// get the context
	auto & context = flecsi::execution::context_t::instance();
	auto rank = context.color();

	constexpr auto num_dims = mesh_t::num_dimensions;
	const auto & vert_lid_to_gid = context.index_map( mesh_t::index_spaces_t::vertices );
	const auto & cell_lid_to_gid = context.index_map( mesh_t::index_spaces_t::cells );

	clog(info) << "DUMP SOLUTION ON RANK " << rank << std::endl;

	// figure out this ranks file name
	auto name_and_ext = ristra::utils::split_extension( filename.str() );
	auto output_filename =
		name_and_ext.first + "_rank" + apps::common::zero_padded(rank) +
		+ "." + name_and_ext.second;

	// dump to file
	if (rank==0)
		std::cout << "Dumping solution to: " << output_filename << std::endl;
	std::ofstream file( output_filename );

	// Dump cell centered quantities
	file.precision(14);
	file.setf( std::ios::scientific );
	file << "# Solution time: " << time << std::endl;
	file << "# Number iterations: " << iteration << std::endl;

	file << "# BEGIN CELLS" << std::endl;
	file << "# Total number: " << mesh.num_cells() << std::endl;
	file << "# local_id global_id ";
	for ( int dim=0; dim<num_dims; ++dim ) file << "centroid(" << dim << ") ";
	file << "density internal_energy pressure ";
	for ( int dim=0; dim<num_dims; ++dim ) file << "velocity(" << dim << ") ";
	file << std::endl;

	for ( auto c : mesh.cells() ) {

		file << c.id() << " " << cell_lid_to_gid.at(c.id()) << " ";
		for ( auto x : c->centroid() ) file << x << " ";
		file << d(c) << " " << e(c) << " " << p(c) << " ";
		for ( auto x : v(c) ) file << x << " ";
		file << std::endl;

	}

	file << "# END CELLS" << std::endl;

	// Dump vertex quantities
	file << "# BEGIN VERTICES" << std::endl;
	file << "# Total number: " << mesh.num_vertices() << std::endl;
	file << "# local_id global_id ";
	for ( int dim=0; dim<num_dims; ++dim ) file << "coordinate(" << dim << ") ";
	file << std::endl;

	for ( auto v : mesh.vertices() ) {

		file << v.id() << " " << vert_lid_to_gid.at(v.id()) << " ";
		for ( auto x : v->coordinates() ) file << x << " ";
		file << std::endl;

	}

	file << "# END VERTICES" << std::endl;

	// close file
	file.close();
}



////////////////////////////////////////////////////////////////////////////////
// TASK REGISTRATION
////////////////////////////////////////////////////////////////////////////////

flecsi_register_task(update_geometry, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(initial_conditions, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(initial_conditions_from_file, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(evaluate_time_step, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(evaluate_gradients, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(evaluate_fluxes, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(apply_update, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(output, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(print, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(dump, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(init_soln_time, apps::hydro, loc, index|flecsi::leaf);
flecsi_register_task(print_soln_time, apps::hydro, loc, index|flecsi::leaf);

} // namespace hydro
} // namespace apps
