/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief A dummy part of the ale library used a template for developpers.
////////////////////////////////////////////////////////////////////////////////


// This is a new way of doing include guards.  You need this in every
// include file to make sure it doesnt get included more than once.
#pragma once 

// Be carefull about including within include files.  Because everone
// that includes this file will also get everything it "#includes"

// system includes could go here
#include <array>

// user includes could also go here ( note the path layout )
//#include "dummy/dummy.h"

// everything in src falls within the ale namespace
namespace ale {

// each subfolder has its own namespace so that we can group functionality
// and easily tell where it "lives"
namespace dummy {

////////////////////////////////////////////////////////////////////////////////
//! \brief function prototypes for functions defined in dummy.cc
////////////////////////////////////////////////////////////////////////////////
int foo();
int bar();


////////////////////////////////////////////////////////////////////////////////
//! \brief some sample class
////////////////////////////////////////////////////////////////////////////////
class foo_t
{
  int i;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief another sample class
////////////////////////////////////////////////////////////////////////////////
class bar_t
{
  int j;
};


} // dummy namespace
} // ale namespace

