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
 * \file ideal_gas.cc
 * 
 * \brief Tests related to the ideal gas.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <cinchtest.h>
#include <iostream>
#include <utility>
#include <vector>

// user includes
#include "ale/eos/ideal_gas.h"
#include "ale/eos/eos_utils.h"


// explicitly use some stuff
using std::cout;
using std::endl;
using std::vector;
using ale::common::real_t;

using namespace ale::eos;

    
///////////////////////////////////////////////////////////////////////////////
//! \brief Test of the ideal gas state class
///////////////////////////////////////////////////////////////////////////////
TEST(eos, ideal_gas) {

  
  ideal_gas_t eos;

  vector<real_t> d{1.0, 2.0};
  vector<real_t> p{1.0, 2.0};

  constexpr size_t n = 2;
  
  constexpr auto rt_two = std::sqrt( 2.0 );

  vector<real_t> e(n), ss(n), t(n);

  compute_internal_energy(eos, d, p, e );
  compute_pressure( eos, d, e, p );
  compute_sound_speed( eos, d, e, ss );
  compute_temperature( eos, d, e, t );

  for ( size_t i = 0; i<n; i++ ) {
    std::cout <<   "d=" <<  d[i] 
              << ", p=" <<  p[i] 
              << ", e=" <<  e[i] 
              << ", s=" << ss[i]
              << ", t=" <<  t[i] << std::endl;
    ASSERT_NEAR( i+1, d[i], TEST_TOLERANCE );
    ASSERT_NEAR( i+1, p[i], TEST_TOLERANCE );
    ASSERT_NEAR( 1.0, e[i], TEST_TOLERANCE );
    ASSERT_NEAR( rt_two, ss[i], TEST_TOLERANCE );
    ASSERT_NEAR( 1.0, t[i], TEST_TOLERANCE );
  }

    
} // TEST_F



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
