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
 * \file tasks.cc
 * 
 * \brief Tests related to generating automatic task functions.
 *
 ******************************************************************************/

// system includes
#include<cinchtest.h>
#include<iostream>
#include<vector>

// user includes
#include "../../utils/tasks.h"

// explicitly use some stuff
using std::cout;
using std::endl;
using std::vector;

using real_t   = double;
using vector_t = vector<real_t>;

using namespace ale::utils;

///////////////////////////////////////////////////////////////////////////////
//! \brief Test the automatic task generator.
///////////////////////////////////////////////////////////////////////////////
TEST(tasks, simple) {

  // create the task
  auto simple_function = [](auto && a, auto && b, auto && c) 
    { c = a + b; };

  vector_t a = { 1.0 }, b = { 2.0 };
  vector_t c = { 0.0 };
  simple_task( simple_function, a, b, c );
  ASSERT_EQ( c[0], 3.0 );
  

} // TEST

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
