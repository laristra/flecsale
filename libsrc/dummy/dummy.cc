/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief A dummy part of the ale library used a template for developpers.
////////////////////////////////////////////////////////////////////////////////

// system includes go here
#include <vector>

// user includes go here ( note the path layout )
#include "dummy/dummy.h"

// everything in src falls within the ale namespace
namespace ale {

// each subfolder has its own namespace so that we can group functionality
// and easily tell where it "lives"
namespace dummy {

////////////////////////////////////////////////////////////////////////////////
//! \brief some sample function
////////////////////////////////////////////////////////////////////////////////
int foo()
{
  foo_t a;
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! \brief another sample function function
////////////////////////////////////////////////////////////////////////////////
int bar()
{
  bar_t b;
  return 0;
}


} // dummy namespace
} // ale namespace
