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
 * \file dummy.cc
 * 
 * \brief Tests a part of the dummy code.
 *
 * This creates an executable test_dummy that runs the tests.    Tests are created
 * with the "TEST(group_name, test_name)" macro.  Executing
 *    ./test_dummy
 * runs all the tests, and calling
 *    ./test_dummy --gtest_list_tests
 * will list all the available tests.  You can execute specific tests by calling
 *    ./test_dummy --gtest_filter=<PATTERN>
 * where <PATTERN> is a typicall shell regex patern.  For example,
 *    ./test_dummy --gtest_filter=dummy.*
 * means run all tests in the dummy group.
 *
 ******************************************************************************/

//
// system includes
//

//! this include is needed for the TEST framework
#include <cinchtest.h>

//! this include is needed for std::cout and std::endl
#include <iostream>

//
// user includes
//

//! this include is needed to get the dummy functionality so we can test it
#include "ale/dummy/dummy.h"


////////////////////////////////////////////////////////////////////////////////
//! \brief This is the first test.
//!
//! The first argument is the test group name, and the second argument
//! is the specific test name.  These names are arbitrary and up to
//! you.  So this test will be called dummy.first_test.  This lets us
//! filter tests --gtest_filter=dummy.* means run all tests in the
//! dummy group.
////////////////////////////////////////////////////////////////////////////////
TEST(dummy, first_test) 
{
  std::cout << "dummy.first_test : HELLO" << std::endl;
  ASSERT_EQ( 0, ale::dummy::foo() );
}


////////////////////////////////////////////////////////////////////////////////
//! \brief Another test
////////////////////////////////////////////////////////////////////////////////
TEST(dummy, second_test) 
{
  std::cout << "dummy.second_test : HELLO" << std::endl;
  ASSERT_FALSE( ale::dummy::bar() );
}

