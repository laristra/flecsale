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
 * \file euler_eqns.cc
 * 
 * \brief Tests related to the euler equations.
 *
 ******************************************************************************/

// system includes
#include <cinchtest.h>
#include <iostream>

// user includes
#include "ale/eqns/euler_eqns.h"
#include "ale/eos/ideal_gas.h"


// explicitly use some stuff
using std::cout;
using std::endl;
using std::vector;
using ale::common::real_t;

using namespace std::placeholders;

using namespace ale::eqns;
using namespace ale::eos;
    
///////////////////////////////////////////////////////////////////////////////
//! \brief Test of the ideal gas state class
///////////////////////////////////////////////////////////////////////////////
TEST(eqns, euler) {

  using eqns_t = euler_eqns_t<3>;
  using eos_t = ideal_gas_t;

  eos_t eos;
  euler_eqns_t<3> eqns;  

  auto get_pressure         = std::bind( &eos_t::compute_pressure,        std::cref(eos), _1, _2 );
  auto get_internal_energy  = std::bind( &eos_t::compute_internal_energy, std::cref(eos), _1, _2 );

  eqns_t::primitive_state_t w{
    1.0, 
    1.0, 0.5,.75, 
    2.0 };
  eqns_t::primitive_state_t w_new;
  eqns_t::conserved_state_t u;

  eqns_t::primitive_to_conserved( w, u, get_pressure );
  eqns_t::conserved_to_primitive( u, w_new, get_internal_energy );
    
} // TEST_F



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
