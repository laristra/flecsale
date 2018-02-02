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
#include <ristra/math/constants.h>
#include <ristra/math/general.h>

namespace apps {
namespace hydro {

// type aliases confined to this translation unit
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
const real_t vol = ristra::math::pi * ristra::math::cube(delta_r) / 6.;

//=============================================================================
// Now set the inputs.
//=============================================================================

// the case prefix
string inputs_t::prefix = "shock_box_3d";
string inputs_t::postfix = "dat";

// output frequency
size_t inputs_t::output_freq = 100;

// the CFL and final solution time
real_t inputs_t::CFL = 1.0/3.0;
real_t inputs_t::final_time = 1.0;
size_t inputs_t::max_steps = 1e6;

// the equation of state
eos_t inputs_t::eos = 
  flecsale::eos::ideal_gas_t<real_t>( 
    /* gamma */ gamma, /* cv */ 1.0 
  ); 

// this is a lambda function to set the initial conditions
inputs_t::ics_function_t inputs_t::ics = 
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

} // namespace
} // namespace
