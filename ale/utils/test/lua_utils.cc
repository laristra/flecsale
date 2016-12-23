/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Tests related to embedded lua.
////////////////////////////////////////////////////////////////////////////////

#if HAVE_LUA

// user includes
#include "ale/common/types.h"
#include "ale/utils/lua_utils.h"

// system includes
#include<cinchtest.h>
#include<iostream>


// explicitly use some stuff
using namespace ale::utils;
using ale::common::test_tolerance;

///////////////////////////////////////////////////////////////////////////////
//! \brief Test the simple use of embedded lua.
///////////////////////////////////////////////////////////////////////////////
TEST(lua_utils, simple) 
{

  // setup the python interpreter
  auto state = lua();

  // load the test file
  state.run_string(
    "print(\"Hello World\")\n"
  );

} // TEST


///////////////////////////////////////////////////////////////////////////////
//! \brief Test the simple use of embedded python.
///////////////////////////////////////////////////////////////////////////////
TEST(lua_utils, embedded) 
{

  // setup the python interpreter
  auto state = lua();

  // load the test file
  state.loadfile( "lua_test.lua" );

  // run a simple function and check the result
  ASSERT_EQ( 3, state.call_function<int>( "sum", 1, 2 ) );
  ASSERT_NEAR( 3., state.call_function<double>( "sum", 1, 2 ), test_tolerance );

  // try with different arguments
  ASSERT_EQ( 3, state.call_function<int>( "sum", 1., 2. ) );
  ASSERT_NEAR( 3., state.call_function<double>( "sum", 1., 2. ), test_tolerance );

} // TEST

#endif // HAVE_LUA
