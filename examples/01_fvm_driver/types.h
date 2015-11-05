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
 * \file types.h
 * 
 * \brief Main types for the example.
 *
 ******************************************************************************/

// flexi includes
#include <flexi/specializations/burton/burton.h>

// select the mesh
using mesh_t = flexi::burton_mesh_t;

// extract some types from the mesh
using vertex_t = flexi::mesh_t::vertex_t;
using   real_t = flexi::mesh_t::  real_t;
using vector_t = flexi::mesh_t::vector_t;


// setup an accessor
using accessor_t =
  flexi::burton_mesh_traits_t::mesh_state_t::accessor_t<real_t>;

//! use some things
using flexi::persistent;
