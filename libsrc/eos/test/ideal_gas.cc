/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Tests related to the ideal gas.
///
////////////////////////////////////////////////////////////////////////////////

// system includes
#include <cinchtest.h>
#include <iostream>
#include <vector>

// user includes
#include "common/types.h"
#include "eos/ideal_gas.h"
#include "utils/tasks.h"


// explicitly use some stuff
using std::cout;
using std::endl;
using std::vector;

using namespace ale;
using namespace ale::eos;
using namespace ale::utils;

using real_t = common::real_t;
using vector_t = vector<real_t>;

using common::test_tolerance;

    
///////////////////////////////////////////////////////////////////////////////
//! \brief Test of the ideal gas state class
///////////////////////////////////////////////////////////////////////////////
TEST(eos, ideal_gas) {


  ideal_gas_t<real_t> eos;

  vector<real_t> d{1.0, 2.0};
  vector<real_t> p{1.0, 2.0};

  constexpr size_t n = 2;
  
  const auto rt_two = std::sqrt( 2.0 );

  vector<real_t> e(n), ss(n), t(n);

  // create the task
  auto update = [&eos](auto && d, auto && p, auto && e, auto && ss, auto && t) 
  { 
    e  = eos.compute_internal_energy_dp( d, p );
    p  = eos.compute_pressure_de( d, e );
    ss = eos.compute_sound_speed_de( d, e );
    t  = eos.compute_temperature_de( d, e );    
  };

  simple_task( update, d, p, e, ss, t );

  // now check
  for ( int i = 0; i<n; i++ ) {
    std::cout <<   "d=" <<  d[i] 
              << ", p=" <<  p[i] 
              << ", e=" <<  e[i] 
              << ", s=" << ss[i]
              << ", t=" <<  t[i] << std::endl;
    ASSERT_NEAR( i+1, d[i], test_tolerance ) << "Density test failed";
    ASSERT_NEAR( i+1, p[i], test_tolerance ) << "Pressure test failed";
    ASSERT_NEAR( 1.0, e[i], test_tolerance ) << "Energy test failed";
    ASSERT_NEAR( rt_two, ss[i], test_tolerance ) << "Sound speed test failed";
    ASSERT_NEAR( 1.0, t[i], test_tolerance ) << "Temperature test failed";
  }

    
} // TEST_F


