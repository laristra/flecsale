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
 * \file mesh.h
 * 
 * \brief Driver functions related to the mesh.
 *
 ******************************************************************************/

#pragma once

// local includes
#include "types.h"

// some mesh constants
constexpr size_t height =  5; //!< number of cells in x
constexpr size_t width  = 11; //!< number of cells in y


////////////////////////////////////////////////////////////////////////////////
//! \brief The main driver for initializing the mesh
//!
//! \param [in,out] mesh the mesh object
//! \param [in] argv an array of command line arguments
////////////////////////////////////////////////////////////////////////////////
int32_t mesh_init(mesh_t & mesh) {

  // initial setup
  auto num_vertices = (height+1)*(width+1);
  mesh.init_parameters( num_vertices );
  
  // create the vertices
  std::vector<vertex_t*> vs;
  for(size_t j = 0; j < height + 1; ++j){
    for(size_t i = 0; i < width + 1; ++i){
      auto v = mesh.create_vertex({real_t(i), real_t(j)});
      v->set_rank(1);
      vs.push_back(v);
    }
  }
  
  // now create the cells
  auto width1 = width + 1;
  for(size_t j = 0; j < height; ++j){
    for(size_t i = 0; i < width; ++i){
      // go over vertices counter clockwise to define cell
      auto c = mesh.create_cell({vs[i + j * width1],
        vs[i + 1 + j * width1],
        vs[i + 1 + (j + 1) * width1],
        vs[i + (j + 1) * width1]});
    } // for
  } // for

  // this sets up everthing else
  mesh.init();

  return 0;

} // driver


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
