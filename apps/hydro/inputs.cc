/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This is where the compile time inputs are defined.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "inputs.h"

namespace apps {
namespace hydro {

// the case prefix
std::string inputs_t::prefix = "sedov_3d";
std::string inputs_t::postfix = "dat";

// output frequency
inputs_t::size_t inputs_t::output_freq = 100;

//! \brief the CFL and final solution time
inputs_t::real_t inputs_t::CFL = 1.0/3.0;
inputs_t::real_t inputs_t::final_time = 0.2;
inputs_t::size_t inputs_t::max_steps = 1e6;

// this is a lambda function to set the initial conditions
inputs_t::ics_function_t inputs_t::ics = []( const inputs_t::vector_t & x )
{
  real_t d, p;
  vector_t v(0);
  if ( x[0] < 0 && x[1] < 0 && x[2] < 0 ) {
    d = 0.125;
    p = 0.1;
  }
  else {
    d = 1.0;
    p = 1.0;
  }    
  return std::make_tuple( d, v, p );
};

// This function builds and returns a mesh
inputs_t::mesh_function_t inputs_t::make_mesh = [](void)
{ 
  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 10;
  constexpr size_t num_cells_z = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 1.0;
  constexpr real_t length_z = 1.0;

  // this is the mesh object
  auto mesh = ale::mesh::box<mesh_t>( 
    num_cells_x, num_cells_y, num_cells_z, length_x, length_y, length_z 
  );

  return mesh;
};


} // namespace
} // namespac
