/*~--------------------------------------------------------------------------~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
/*!
 *
 * \file fvm_example.h
 * 
 * \brief An example driver for a simple finite-volume solver.
 *
 ******************************************************************************/

#pragma once

// system includes
#include <iostream>

// flexi includes
#include <flecsi/execution/task.h>

// local includes
#include "mesh.h"



////////////////////////////////////////////////////////////////////////////////
//! \brief The main driver for the finite-volume solver
//!
//! \param [in] argc the number of command line arguments
//! \param [in] argv an array of command line arguments
////////////////////////////////////////////////////////////////////////////////
int32_t driver(int argc, char ** argv) {

  mesh_t mesh;
  //equation_t eqn;
  //discrete_t<equation_t> disc;

  // execute a task to initialize the mesh
  execute(mesh_init, mesh);

  // register several state variables
  register_state(mesh,  "density", cells,   real_t, persistent);
  register_state(mesh, "pressure", cells,   real_t, persistent);
  register_state(mesh, "velocity", cells, vector_t, persistent);

  // register several state variables  
  //execute(init_solver, mesh);

  // initialize the solution
  //execute(apply_ics, mesh);

  // intial write
  //execute(write, mesh);

  // execute solver
  //execute(solve, mesh);

  // final write
  //execute(write, mesh);

  return 0;

}

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
