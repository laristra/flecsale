/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// 
/// \file
///
/// \brief This header defines some usefull types.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <cstddef> // for size_t
#include <cstdint> // for types

namespace flecsale {
namespace common {

//! function return values
enum return_t { 
  RETURN_SUCCESS = 0, //! command succeeded 
  RETURN_FAILURE = 1  //! command failed
};

//! real precision type
#ifdef FleCSALE_DOUBLE_PRECISION
using real_t = double;
#else
using real_t = float;
#endif

//! type of integer ids to use
#ifdef FleCSALE_USE_64BIT_IDS
using index_t = uint64_t;
#else
using index_t = uint32_t;
#endif

//! type of integer data to use
#ifdef FleCSALE_DOUBLE_PRECISION
using integer_t = int64_t;
#else
using integer_t = int32_t;
#endif

//! a size type
using std::size_t;

//! a counter type
using counter_t = long long;

} // namespace
} // namespace


