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
 * \file factory.h
 * 
 * \brief Some functionality for creating meshes.
 *
 ******************************************************************************/

#pragma once

namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \brief Create a box mesh
//!
//! \param [in] num_cells_x,num_cells_y  the number of cells in the x and y dir
//! \param [in] min_x,min_y              the min coordinate in the x and y dir
//! \param [in] max_x,max_y              the max coordinate in the x and y dir
//! \return a new mesh object
////////////////////////////////////////////////////////////////////////////////
template< typename T >
T box( typename T::size_t num_cells_x, typename T::size_t num_cells_y, 
       typename T::real_t min_x,       typename T::real_t min_y,
       typename T::real_t max_x,       typename T::real_t max_y ) 
{

  T mesh;

  // the grid dimensions
  auto length_x = max_x - min_x;
  auto length_y = max_y - min_y;


  // reserve storage for the mesh
  auto num_vertex = ( num_cells_x + 1 ) * ( num_cells_y + 1 );
  mesh.init_parameters( num_vertex );
  
  
  // create the individual vertices
  using vertex_t = typename T::vertex_t;
  std::vector<vertex_t*> vs;
  
  auto delta_x = length_x / num_cells_x;
  auto delta_y = length_y / num_cells_y;

  auto num_vert_x = num_cells_x + 1;
  auto num_vert_y = num_cells_y + 1;

  for(size_t j = 0; j < num_vert_y; ++j) {
    auto y = min_y + j*delta_y;
    for(size_t i = 0; i < num_vert_x; ++i) {
      auto x = min_x + i*delta_x;
      auto v = mesh.create_vertex( {x, y} );
      vs.push_back(v);
    }
    
  }
  
  // define each cell
  auto index = [=](auto i, auto j) { return i + num_vert_x*j; };
  
  for(size_t j = 0; j < num_cells_y; ++j)
    for(size_t i = 0; i < num_cells_x; ++i) {
      auto c = 
        mesh.create_cell({
              vs[ index(i  , j  ) ],
              vs[ index(i+1, j  ) ],
              vs[ index(i+1, j+1) ],
              vs[ index(i  , j+1) ]});
    }
  
  
  // now finalize the mesh setup
  mesh.init();

  return mesh;
}



////////////////////////////////////////////////////////////////////////////////
//! \brief Create a box mesh
//!
//! \param [in] num_cells_x,num_cells_y  the number of cells in the x and y dir
//! \param [in] min_x,min_y              the min coordinate in the x and y dir
//! \param [in] max_x,max_y              the max coordinate in the x and y dir
//! \return a new mesh object
////////////////////////////////////////////////////////////////////////////////
template< typename T >
T box( typename T::size_t num_cells_x, typename T::size_t num_cells_y, 
       typename T::real_t length_x,    typename T::real_t length_y ) 
{
  
  auto x1 =  length_x / 2;
  auto y1 =  length_y / 2;
  auto x0 = - x1;
  auto y0 = - y1;

  return box<T>( num_cells_y, num_cells_x, x0, y0, x1, y1 );
}

} // namespace mesh
} // namespace ale
