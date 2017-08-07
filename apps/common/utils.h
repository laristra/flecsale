/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some commonly used utilities.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <sstream>

namespace apps {
namespace common {

///////////////////////////////////////////////////////////////////////////////
//! \brief Tack on an iteration number to a string
///////////////////////////////////////////////////////////////////////////////
static auto zero_padded( 
  std::size_t n, std::size_t padding = 6 
)
{
  std::stringstream ss;
  ss << std::setw( padding ) << std::setfill( '0' ) << n;
  return ss.str();;
}

} // namespace
} // namespace
