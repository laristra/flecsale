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
int32_t initial_conditions_task( 
  mesh::burton_mesh_3d_t & mesh, inputs_t::ics_function_t ics 
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
int32_t update_state_from_pressure_task( mesh::burton_mesh_3d_t & mesh ) 
{
	return update_state_from_pressure( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task for updating the state from energy.
//!
//! Updates the state from density and energy and computes the new pressure.
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t update_state_from_energy_task( mesh::burton_mesh_3d_t & mesh ) 
{
	return update_state_from_energy( mesh );
}


////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute the time step size
//!
//! \param [in,out] mesh  the mesh object
//! \param [in,out] limit_string  a string describing the limiting time step
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t evaluate_time_step_task( 
	mesh::burton_mesh_3d_t & mesh, std::string & limit_string 
) {
  return evaluate_time_step( mesh, limit_string );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t estimate_nodal_state_task( mesh::burton_mesh_3d_t & mesh ) 
{
  return estimate_nodal_state( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t evaluate_corner_coef_task( mesh::burton_mesh_3d_t & mesh ) 
{
  return evaluate_corner_coef( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to compute nodal quantities
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t evaluate_nodal_state_task( 
  mesh::burton_mesh_3d_t & mesh, 
  const boundary_map_t<mesh::burton_mesh_3d_t::num_dimensions> & boundary_map
) {
  return evaluate_nodal_state( mesh, boundary_map );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to evaluate residuals
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t evaluate_forces_task( mesh::burton_mesh_3d_t & mesh ) 
{
  return evaluate_forces( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to update the solution
//!
//! \param [in,out] mesh the mesh object
//!   \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t apply_update_task( 
  mesh::burton_mesh_3d_t & mesh, real_t coef, bool first_time 
) {
  return apply_update( mesh, coef, first_time );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to move the mesh
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t move_mesh_task( mesh::burton_mesh_3d_t & mesh, real_t coef ) 
{
  return move_mesh( mesh, coef );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t save_coordinates_task( mesh::burton_mesh_3d_t & mesh ) 
{
  return save_coordinates( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t restore_coordinates_task( mesh::burton_mesh_3d_t & mesh ) 
{
  return restore_coordinates( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to save the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t save_solution_task( mesh::burton_mesh_3d_t & mesh ) 
{
  return save_solution( mesh );
}

////////////////////////////////////////////////////////////////////////////////
//! \brief The main task to restore the coordinates
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
int32_t restore_solution_task( mesh::burton_mesh_3d_t & mesh ) 
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
flecsi_register_task(estimate_nodal_state_task, loc, single);
flecsi_register_task(evaluate_corner_coef_task, loc, single);
flecsi_register_task(evaluate_nodal_state_task, loc, single);
flecsi_register_task(evaluate_forces_task, loc, single);
flecsi_register_task(apply_update_task, loc, single);
flecsi_register_task(move_mesh_task, loc, single);
flecsi_register_task(save_coordinates_task, loc, single);
flecsi_register_task(restore_coordinates_task, loc, single);
flecsi_register_task(save_solution_task, loc, single);
flecsi_register_task(restore_solution_task, loc, single);

} // namespace
} // namespace
