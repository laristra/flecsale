/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "inputs.h"
#include "tasks.h"
#include "../driver.h"


namespace flecsi {
namespace execution {


///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
void driver(int argc, char** argv) 
{
  using namespace apps::hydro;
  apps::hydro::driver<inputs_t>( argc, argv );
}

} // namespace
} // namespace
