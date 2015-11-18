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
  eqns_t eqns;  

  using   real_t = eqns_t::real_t;
  using vector_t = eqns_t::vector_t;

  auto get_pressure         = std::bind( &eos_t::compute_pressure,        std::cref(eos), _1, _2 );
  auto get_internal_energy  = std::bind( &eos_t::compute_internal_energy, std::cref(eos), _1, _2 );

  eqns_t::primitive_state_t w{
    real_t{1.0}, 
    vector_t{1.0, 0.5,.75}, 
    real_t{2.0} };

  eqns_t::conserved_state_t u     = w.to_conserved( get_pressure );
  eqns_t::primitive_state_t w_new = u.to_primitive( get_internal_energy );

  eqns_t::conserved_state_t u_cpy(u);

  cout << w     << endl;
  cout << u     << endl;
  cout << w_new << endl;

  ASSERT_TRUE( w == w_new );
  ASSERT_TRUE( u == u_cpy );

    
} // TEST_F



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
