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

namespace apps {
namespace hydro {

// type aliases confined to this translation unit
using size_t = inputs_t::size_t;
using real_t = inputs_t::real_t;
using vector_t = inputs_t::vector_t;
using string = std::string;
using eos_t = inputs_t::eos_t;


// the case prefix
string inputs_t::prefix = "shock_box_2d";
string inputs_t::postfix = "dat";

// output frequency
size_t inputs_t::output_freq = 1e6;

// the CFL and final solution time
real_t inputs_t::CFL = 1.0/2.0;
real_t inputs_t::final_time = 0.2;
size_t inputs_t::max_steps = 20;

// the equation of state
eos_t inputs_t::eos = 
  flecsale::eos::ideal_gas_t<real_t>( 
    /* gamma */ 1.4, /* cv */ 1.0 
  ); 

// this is a lambda function to set the initial conditions
inputs_t::ics_function_t inputs_t::ics = 
  []( const vector_t & x, const real_t & )
  {
    real_t d, p;
    vector_t v(0);
    if ( x[0] < 0 && x[1] < 0 ) {
      d = 0.125;
      p = 0.1;
    }
    else {
      d = 1.0;
      p = 1.0;
    }    
    return std::make_tuple( d, v, p );
  };

} // namespace
} // namespace
