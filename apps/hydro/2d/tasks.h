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
#include "../types.h"
#include "../tasks.h"

#include <flecsi/execution/context.h>
#include <flecsi/execution/execution.h>

namespace apps {
namespace hydro {

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for setting initial conditions
//!
//! \param [in,out] mesh the mesh object
//! \param [in]     ics  the initial conditions to set
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int initial_conditions_task( 
  mesh::burton_mesh_2d_t & mesh, inputs_t::ics_function_t ics 
) {
  return initial_conditions( mesh, ics );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for updating the state using pressure.
//!
//! Updates the state from density and pressure and computes the new energy.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int update_state_from_pressure_task( 
  const mesh::burton_mesh_2d_t & mesh, const eos_t * eos
) {
	return update_state_from_pressure( mesh, eos );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for updating the state from energy.
//!
//! Updates the state from density and energy and computes the new pressure.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int update_state_from_energy_task( 
  mesh::burton_mesh_2d_t & mesh, const eos_t * eos
) {
	return update_state_from_energy( mesh, eos );
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int evaluate_time_step_task( mesh::burton_mesh_2d_t & mesh ) 
{
  using eqns_t = eqns_t<mesh::burton_mesh_2d_t::num_dimensions>;
  return evaluate_time_step<eqns_t>( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate fluxes at each face.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int evaluate_fluxes_task( mesh::burton_mesh_2d_t & mesh ) 
{
  return evaluate_fluxes( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution in each cell.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
solution_error_t apply_update_task( 
  mesh::burton_mesh_2d_t & mesh, real_t tolerance, bool first_time
) {
  return apply_update( mesh, tolerance, first_time );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int save_solution_task( mesh::burton_mesh_2d_t & mesh ) 
{
  return save_solution( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int restore_solution_task( mesh::burton_mesh_2d_t & mesh ) 
{
  return restore_solution( mesh );
}

////////////////////////////////////////////////////////////////////////////////
// TASK REGISTRATION
////////////////////////////////////////////////////////////////////////////////

flecsi_register_task(initial_conditions_task, loc, single);
flecsi_register_task(update_state_from_pressure_task, loc, single);
flecsi_register_task(update_state_from_energy_task, loc, single);
flecsi_register_task(evaluate_time_step_task, loc, single);
flecsi_register_task(evaluate_fluxes_task, loc, single);
flecsi_register_task(apply_update_task, loc, single);
flecsi_register_task(save_solution_task, loc, single);
flecsi_register_task(restore_solution_task, loc, single);

} // namespace
} // namespace
