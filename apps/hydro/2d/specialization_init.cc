/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "types.h"
#include "../specialization_init.h"

namespace flecsi {
namespace execution {

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver.
///////////////////////////////////////////////////////////////////////////////
void specialization_tlt_init(int argc, char** argv) 
{
  apps::hydro::specialization_tlt_init( argc, argv );
}

} // namespace
} // namespace
