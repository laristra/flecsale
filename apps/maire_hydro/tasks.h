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
#include "globals.h"
#include "types.h"

#include <flecsale/io/io_exodus.h>
#include <flecsale/linalg/qr.h>
#include <ristra/utils/algorithm.h>
#include <ristra/utils/array_view.h>
#include <ristra/utils/filter_iterator.h>
#include <ristra/utils/string_utils.h>
#include <flecsi/execution/context.h>
#include <flecsi/execution/execution.h>

// system includes
#include <iomanip>

namespace apps {
namespace hydro {


////////////////////////////////////////////////////////////////////////////////
//! \brief Install a boundary on the mesh
////////////////////////////////////////////////////////////////////////////////
void install_boundary( 
  client_handle_r__<mesh_t>  mesh,
  real_t soln_time,
  tag_t bc_key,
  boundary_condition_t * bc_type,
  inputs_t::bcs_function_t bc_function
) {

  mesh.install_boundary( 
      [=](auto f) 
      { 
        if ( f->is_boundary() ) {
          const auto & fx = f->midpoint();
          return ( bc_function(fx, soln_time) );
        }
        return false;
      },
      bc_key
    );
    globals::boundaries.emplace( bc_key, bc_type );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Check if the mesh is correct
//!
//! \param [in] mesh the mesh object
////////////////////////////////////////////////////////////////////////////////
void validate_mesh( 
  client_handle_r__<mesh_t> mesh
) {
  
  mesh.is_valid();

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void initial_conditions( 
  client_handle_r__<mesh_t>  mesh,
  inputs_t::ics_function_t ics, 
  eos_t eos,
  real_t soln_time,
  dense_handle_w__<real_t> V,
  dense_handle_w__<real_t> M,
  dense_handle_w__<vector_t> v,
  dense_handle_w__<real_t> p,
  dense_handle_w__<real_t> d,
  dense_handle_w__<real_t> e,
  dense_handle_w__<real_t> T,
  dense_handle_w__<real_t> a
) {

  // This doesn't work with lua input
  // #pragma omp parallel for
  for ( auto c : mesh.cells( flecsi::owned ) ) {
    // now copy the state to flexi
    real_t den;
    std::tie( den, v(c), p(c) ) = ics( c->centroid(), soln_time );
    // set mass and volume now
    const auto & cell_vol = c->volume();
    M(c) = den*cell_vol;
    V(c) = cell_vol;
    // now update the rest of the state
    eqns_t::update_state_from_pressure( 
      pack( c, V, M, v, p, d, e, T, a ),
      eos
    );
  }

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void update_state_from_energy( 
  client_handle_r__<mesh_t>  mesh,
  eos_t eos,
  dense_handle_r__<real_t> V,
  dense_handle_r__<real_t> M,
  dense_handle_r__<vector_t> v,
  dense_handle_w__<real_t> p,
  dense_handle_r__<real_t> d,
  dense_handle_r__<real_t> e,
  dense_handle_w__<real_t> T,
  dense_handle_w__<real_t> a
) {

  auto cs = mesh.cells( flecsi::owned );
  auto num_cells = cs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_cells; ++i ) {
    auto c = cs[i];
    auto u = pack(c, V, M, v, p, d, e, T, a);
    eqns_t::update_state_from_energy( u, eos );
  }

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size
//!
//! \param [in,out] mesh  the mesh object
//! \param [in,out] limit_string  a string describing the limiting time step
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
double evaluate_time_step(
  client_handle_r__<mesh_t> mesh,
  time_constants_t cfl, 
	real_t previous_time_step,
	dense_handle_r__<real_t> sound_speed,
	dense_handle_r__<flux_data_t> dudt
) {
 
  // Loop over each cell, computing the minimum time step,
  // which is also the maximum 1/dt
  real_t dt_acc_inv(0);
  real_t dt_vol_inv(0);

  auto cs = mesh.cells( flecsi::owned );
  auto num_cells = cs.size();

  for ( counter_t i=0; i<num_cells; ++i ) {
    auto c = cs[i];

    // compute the inverse of the time scale
    auto dti =  sound_speed(c) / c->min_length();
    // check for the maximum value
    dt_acc_inv = std::max( dti, dt_acc_inv );

    // now check the volume change
    auto dVdt = eqns_t::volumetric_rate_of_change( dudt(c) );
    dti = std::abs(dVdt) / c->volume();
    // check for the maximum value
    dt_vol_inv = std::max( dti, dt_vol_inv );

  } // cell


  assert( dt_acc_inv > 0 && "infinite delta t" );
  assert( dt_vol_inv > 0 && "infinite delta t" );
  
  // get the individual cfls
  auto dt_acc = cfl.accoustic / dt_acc_inv;
  auto dt_vol = cfl.volume / dt_vol_inv;
  auto dt_growth = cfl.growth * previous_time_step;

  // find the minimum one
  auto dts = std::array<real_t,3> { dt_acc, dt_vol, dt_growth };
  auto it = std::min_element( dts.begin(), dts.end() );

	// return the solution
  return *it;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void estimate_nodal_state( 
  client_handle_r__<mesh_t>  mesh,
  dense_handle_r__<vector_t> cell_vel,
  dense_handle_w__<vector_t> vertex_vel // Hack to avoid communication
) {

  using subset_t = mesh_t::subset_t;
  for ( auto v : mesh.vertices(subset_t::overlapping) )
  {
    vertex_vel(v) = 0.;
    const auto & cells = mesh.cells(v);
    for ( auto c : cells ) vertex_vel(v) += cell_vel(c);
    vertex_vel(v) /= cells.size();
  } // vertex

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void evaluate_nodal_state( 
  client_handle_r__<mesh_t>  mesh,
  real_t soln_time,
  dense_handle_r__<real_t> Vc,
  dense_handle_r__<real_t> Mc,
  dense_handle_r__<vector_t> uc,
  dense_handle_r__<real_t> pc,
  dense_handle_r__<real_t> dc,
  dense_handle_r__<real_t> ec,
  dense_handle_r__<real_t> Tc,
  dense_handle_r__<real_t> ac,
  dense_handle_w__<vector_t> un,
  dense_handle_w__<vector_t> npc,
  dense_handle_w__<vector_t> Fpc
) {

  // get the number of dimensions and create a matrix
  constexpr auto num_dims = mesh_t::num_dimensions;

  // use a matrix type
  using matrix_t = ristra::math::matrix<real_t, num_dims, num_dims>;
  
  // the subsets type
  using subset_t = mesh_t::subset_t;

  //----------------------------------------------------------------------------
  // Loop over each vertex
  //----------------------------------------------------------------------------
  for ( auto vt : mesh.vertices( subset_t::overlapping ) ) {

    // create the final matrix the point
    matrix_t Mp(0);
    vector_t rhs(0);

    // get the corners
    auto cnrs = mesh.corners(vt);
    auto num_corners = cnrs.size();

    // create some corner storage
    std::vector< matrix_t > Mpc(num_corners, 0);

    //--------------------------------------------------------------------------
    // build point matrix
    for ( int j=0; j<num_corners; ++j ) {

      // get the corner
      auto cn = cnrs[j];

      // initialize the corner force
      Fpc(cn) = 0;
      npc(cn) = 0;

      // corner attaches to one cell and one point
      auto cl = mesh.cells(cn).front();
      // get the cell state (there is only one)
      auto state = pack(cl, Vc, Mc, uc, pc, dc, ec, Tc, ac);

      // the corner quantities are approximated as cell ones
      const auto & pc = eqns_t::pressure( state );
      const auto & uc = eqns_t::velocity( state );
      const auto & dc = eqns_t::density( state );
      const auto & ac = eqns_t::sound_speed( state );
      // the corner impedance
      auto zc = dc * ac;

      // iterate over the wedges in pairs
      auto ws = mesh.wedges(cn);
      for ( auto w : ws ) 
      {
        // get the first wedge normal
        const auto & n = w->facet_normal();
        const auto & l = w->facet_area();
        // the final matrix
        // Mpc = zc * ( lpc^- npc^-.npc^-  + lpc^+ npc^+.npc^+ );
        ristra::math::outer_product( n, n, Mpc[j], zc*l );
        // compute the pressure coefficient
        for ( int d=0; d<num_dims; ++d ) 
          npc(cn)[d] += l * n[d];
      } // wedges

      // add to the global matrix
      Mp += Mpc[j];
      // compute a portion of the corner force and 
      // add the pressure and velocity contributions to the system
      ax_plus_y( Mpc[j], uc, Fpc(cn) );   
      for ( int d=0; d<num_dims; ++d ) {
        Fpc(cn)[d] += pc * npc(cn)[d];
        rhs[d] += Fpc(cn)[d];
      }

    } // corner

    //--------------------------------------------------------------------------
    // now solve the system for the point velocity

    //---------- boundary point
    if ( vt->is_boundary() ) {

      // this is used to keep track of the symmetry normals
      std::map< tag_t, vector_t > symmetry_normals;

      // get the boundary tags
      const auto & point_tags =  vt->tags();

      // first check if this has a prescribed velocity.  If it does, then nothing to do
      auto vel_bc = std::find_if( 
        point_tags.begin(), point_tags.end(), 
        [&boundaries=globals::boundaries](const auto & id) { 
          return boundaries.at(id)->has_prescribed_velocity();
        } );
      if ( vel_bc != point_tags.end() ) {
        un(vt) = 
          globals::boundaries.at(*vel_bc)->velocity( vt->coordinates(), soln_time );
        continue;
      }

      // otherwise, apply the pressure conditions
      for ( auto w : mesh.wedges(vt) ) 
      {
				// skip internal wedges
				if ( ! w->is_boundary() ) continue;
        // for wedges attached to boundary
				auto f = mesh.faces(w).front();
        auto c = mesh.cells(w).front();
        for ( auto tag : f->tags() ) {
          auto b = globals::boundaries.at( tag );
          // PRESSURE CONDITION
          if ( b->has_prescribed_pressure() ) {
            const auto & n = w->facet_normal();
            const auto & l = w->facet_area();
            const auto & x = w->facet_centroid();
            auto fact = l * b->pressure( x, soln_time );
            for ( int d=0; d<num_dims; ++d )
              rhs[d] -= fact * n[d];
          }
          // SYMMETRY CONDITION
          else if ( b->has_symmetry() ) {
            const auto & n = w->facet_normal();
            const auto & l = w->facet_area();
            if ( symmetry_normals.count(tag) ) {
              auto & tmp = symmetry_normals[tag];
              for ( int d=0; d<num_dims; ++d )
                tmp[d] += l * n[d];
            }
            else {
              auto & tmp = symmetry_normals[tag];
              for ( int d=0; d<num_dims; ++d )              
                tmp[d] = l * n[d];
            }
          } // END CONDITIONS
        } // for each tag
      } // for each wedge


      // now construct the system to solve
      //
      // no additional symmetry constraints
      if ( symmetry_normals.empty() ) {
        un(vt)  = ristra::math::solve( Mp, rhs );
      }
      // add symmetry constraints and grow the system
      else {
        // how many extra constraints ( there are two wedges per face )
        auto num_symmetry = symmetry_normals.size();
        // the matrix size
        auto num_rows = num_dims+num_symmetry;
        // create storage for the new system in a 1d array
        std::vector< real_t > A( num_rows * num_rows, 0 ); // zerod
        std::vector< real_t > b( num_rows ); // zerod
        // create the views
        auto A_view = ristra::utils::make_array_view( A, num_rows, num_rows );
        auto M_view = ristra::utils::make_array_view( Mp.data(), num_dims, num_dims );
        auto b_view = ristra::utils::make_array_view( b );
        // insert the old system into the new one
        for ( int d=0; d<num_dims; ++d )
          b_view[d] = rhs[d];
        for ( int i=0; i<num_dims; i++ ) 
          for ( int j=0; j<num_dims; j++ ) 
            A_view(i,j) = Mp(i,j);
        // insert each constraint
        for ( int i=0; i<num_dims; i++ ) {
          int j = num_dims;
          for ( const auto & n : symmetry_normals )
            A_view( i, j++ ) = n.second[i];          
        }
        int i = num_dims;
        for ( const auto & n : symmetry_normals ) {
          for ( int j=0; j<num_dims; j++ ) 
            A_view( i, j ) = n.second[j];          
          i++;
        }               
        // solve the system
        flecsale::linalg::qr( A_view, b_view );
        // copy the results back
        for ( int d=0; d<num_dims; ++d )
          un(vt)[d] = b_view[d];
      } // end has symmetry

    } // boundary point

    //---------- internal point
    // make sure sum(lpc) = 0
    // assert( abs(np) < eps && "error in norms" );
    // now solve for point velocity
    else {
      
      un(vt) = ristra::math::solve( Mp, rhs );

    } // internal point


    //--------------------------------------------------------------------------
    // Scatter RHS
    for ( int j=0; j<num_corners; ++j ) {

      // get the corner
      auto cn = cnrs[j];
      // corner attaches to one cell and one point
      auto cl = mesh.cells(cn).front();
    
      // now add the vertex component to the force
      matrix_vector( 
        static_cast<real_t>(-1), Mpc[j], un(vt), 
        static_cast<real_t>(1), Fpc(cn)
      );

    }

  } // vertex
  //----------------------------------------------------------------------------

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to sum the forces and comput the cell changes
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void evaluate_residual( 
  client_handle_r__<mesh_t>  mesh,
  dense_handle_r__<vector_t> uv,
  dense_handle_r__<vector_t> npc,
  dense_handle_r__<vector_t> Fpc,
  dense_handle_w__<flux_data_t> dudt // hack so no communication occurs
)
{

  // TASK: loop over each cell and compute the residual

  for ( auto cl : mesh.cells(flecsi::owned) ) {
    
    // Gather corner forces to compute the cell residual

    // local cell residual
    dudt(cl) = 0;

    // compute subcell forces
    for ( auto cn : mesh.corners(cl) ) {
      // corner attaches to one point and zone
      auto pt = mesh.vertices(cn).front();
      // add contribution
      eqns_t::compute_update( uv(pt), Fpc(cn), npc(cn), dudt(cl) );
    }// corners    
    
  } // cell
    
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution
//!
//! \param [in,out] mesh the mesh object
//!   \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void apply_update(
  client_handle_r__<mesh_t>  mesh,
	real_t delta_t,
  dense_handle_r__<flux_data_t> dudt,
  dense_handle_w__<real_t> Vc,
  dense_handle_r__<real_t> Mc,
  dense_handle_w__<vector_t> uc,
  dense_handle_r__<real_t> pc,
  dense_handle_w__<real_t> dc,
  dense_handle_w__<real_t> ec,
  dense_handle_r__<real_t> Tc,
  dense_handle_r__<real_t> ac
) {

  // Using the cell residual, update the state
  for ( auto cl : mesh.cells(flecsi::owned) ) {

    // get the cell state
    auto u = pack(cl, Vc, Mc, uc, pc, dc, ec, Tc, ac);

    // apply the update
    eqns_t::update_state_from_flux( u, dudt(cl), delta_t );
    eqns_t::update_volume( u, cl->volume() );

  } // for

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to move the mesh
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void move_mesh(
	 client_handle_r__<mesh_t> mesh,
	 dense_handle_r__<vector_t> vel,
	 real_t delta_t
) {

  // Update ALL vertices, including ghost so that we dont need to communicate.
	// DEFECT we are modifying the mesh, but its read-only.

  for ( auto vt : mesh.vertices() ) {
    for ( int d=0; d<mesh_t::num_dimensions; ++d )
      vt->coordinates()[d] += delta_t * vel(vt)[d];
  }

	// now update the geometry
	mesh.update_geometry();

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void save_coordinates( 
  client_handle_r__<mesh_t>  mesh,
  dense_handle_r__<vector_t> coord0 // Hack to avoid communication
)
{

  // Loop over vertices
  auto vs = mesh.vertices();
  auto num_verts = vs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_verts; i++ ) {
    auto vt = vs[i];
    coord0(vt) = vt->coordinates();
  }

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void restore_coordinates( 
  client_handle_r__<mesh_t>  mesh,
  dense_handle_r__<vector_t> coord0
)
{

  // Loop over vertices
  auto vs = mesh.vertices();
  auto num_verts = vs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_verts; i++ ) {
    auto vt = vs[i];
    vt->coordinates() = coord0(vt);
  }

}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void save_solution( 
  client_handle_r__<mesh_t>  mesh,
  dense_handle_r__<vector_t> cell_vel,
	dense_handle_r__<real_t> cell_ener,
	dense_handle_r__<vector_t> cell_vel_0, // "r" permissions hack to avoid communication
	dense_handle_r__<real_t> cell_ener_0   // "r" permissions hack to avoid communication
)
{

  // Loop over cells
  auto cs = mesh.cells();
  auto num_cells = cs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_cells; i++ ) {
    auto c = cs[i];
    cell_vel_0(c) = cell_vel(c);
    cell_ener_0(c) = cell_ener(c);
  }

}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
void restore_solution( 
  client_handle_r__<mesh_t>  mesh,
	dense_handle_r__<vector_t> cell_vel_0,
  dense_handle_w__<vector_t> cell_vel,
	dense_handle_r__<real_t> cell_ener_0,
	dense_handle_w__<real_t> cell_ener
)
{

  // Loop over cells
  auto cs = mesh.cells();
  auto num_cells = cs.size();

  #pragma omp parallel for
  for ( counter_t i=0; i<num_cells; i++ ) {
    auto c = cs[i];
    cell_vel(c) = cell_vel_0(c);
    cell_ener(c) = cell_ener_0(c);
  }

}



////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
////////////////////////////////////////////////////////////////////////////////
void output( 
  client_handle_r__<mesh_t> mesh, 
  char_array_t prefix,
	char_array_t postfix,
	size_t iteration,
	real_t time,
  dense_handle_r__<real_t> d,
  dense_handle_r__<vector_t> v,
  dense_handle_r__<real_t> e,
  dense_handle_r__<real_t> p,
  dense_handle_r__<real_t> T,
  dense_handle_r__<real_t> a
) {
  clog(info) << "OUTPUT MESH TASK" << std::endl;
 
  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // figure out this ranks file name
  auto output_filename = 
    prefix.str() + "_rank" + apps::common::zero_padded(rank) +
    "_" + apps::common::zero_padded(iteration) + "." + postfix.str();

  // now outut the mesh
  flecsale::io::io_exodus__<mesh_t>::write(
    output_filename, mesh, iteration, time, &d //, v, e, p, T, a
  );
}

////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
////////////////////////////////////////////////////////////////////////////////
void print( 
  client_handle_r__<mesh_t> mesh,
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
// TASK REGISTRATION
////////////////////////////////////////////////////////////////////////////////

flecsi_register_task(validate_mesh, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(initial_conditions, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(install_boundary, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(estimate_nodal_state, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(evaluate_nodal_state, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(evaluate_residual, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(evaluate_time_step, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(move_mesh, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(apply_update, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(update_state_from_energy, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(save_coordinates, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(restore_coordinates, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(save_solution, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(restore_solution, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(output, apps::hydro, loc, single|flecsi::leaf);
flecsi_register_task(print, apps::hydro, loc, single|flecsi::leaf);

} // namespace hydro
} // namespace apps
