// system includes
#include <cinchtest.h>
#include <iostream>
#include <vector>

// user includes
#include "ale/eos/ideal_gas.h"


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

  eos.compute_internal_energy( d, p, e );
  eos.compute_pressure( d, e, p );
  eos.compute_sound_speed( d, e, ss );
  eos.compute_temperature( d, e, t );

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
