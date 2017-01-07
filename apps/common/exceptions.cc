/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utilities for exception handling.
////////////////////////////////////////////////////////////////////////////////


// system includes
#ifdef _GNU_SOURCE
#  include <fenv.h>
#endif

///////////////////////////////////////////////////////////////////////////////
// enable exceptions
///////////////////////////////////////////////////////////////////////////////
void enable_exceptions(void) {

  // enable exceptions
#if defined(_GNU_SOURCE) && !defined(NDEBUG)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

}
