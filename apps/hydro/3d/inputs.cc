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

#include <flecsale/eos/ideal_gas.h>
#include <flecsale/math/constants.h>
#include <flecsale/math/general.h>

namespace apps {
namespace hydro {

// type aliases confined to this translation unit
using base_t = inputs_t::base_t;
using size_t = inputs_t::size_t;
using real_t = inputs_t::real_t;
using vector_t = inputs_t::vector_t;
using string = std::string;
using eos_t = inputs_t::eos_t;

//=============================================================================
// These constants are not part of the input class, they are globally scoped
// within this translation unit only.
//=============================================================================

// the value of gamma
const real_t gamma = 1.4;

// compute a radial size
const real_t delta_r = .1;

// and a representative volume
// - 1/8-th the volume of a sphere
const real_t vol = math::pi * math::cube(delta_r) / 6.;

//=============================================================================
// Now set the inputs.
//=============================================================================

// the case prefix
template<> string base_t::prefix = "shock_box_3d";
template<> string base_t::postfix = "dat";

// output frequency
template<> size_t base_t::output_freq = 100;

// the CFL and final solution time
template<> real_t base_t::CFL = 1.0/3.0;
template<> real_t base_t::final_time = 1.0;
template<> size_t base_t::max_steps = 1e6;

// the equation of state
template<> eos_t base_t::eos = 
  flecsale::eos::ideal_gas_t<real_t>( 
    /* gamma */ gamma, /* cv */ 1.0 
  ); 

// this is a lambda function to set the initial conditions
template<>
inputs_t::ics_function_t base_t::ics = 
  [g=gamma]( const vector_t & x, const real_t & )
  {
    constexpr real_t e0 = 0.106384;
    real_t d = 1.0;
    vector_t v = 0;
    real_t p = 1.e-6;
    auto r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
    if ( r < delta_r + std::numeric_limits<real_t>::epsilon() )
      p = (g - 1) * d * e0 / vol;
    return std::make_tuple( d, v, p );
  };

// This function builds and returns a mesh
template<>
inputs_t::mesh_function_t base_t::make_mesh = 
  [](const real_t &)
  {
    // the grid dimensions
    constexpr size_t num_cells_x = 10;
    constexpr size_t num_cells_y = 10;
    constexpr size_t num_cells_z = 10;
    constexpr real_t length_x = 1.0;
    constexpr real_t length_y = 1.0;
    constexpr real_t length_z = 1.0;
  
    // this is the mesh object
    return flecsale::mesh::box<mesh_t>( 
      num_cells_x, num_cells_y, num_cells_z, length_x, length_y, length_z
    );
  };

} // namespace
} // namespace
