/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines some 2d-specific types.
///////////////////////////////////////////////////////////////////////////////
#pragma once 

// user includes
#include "../types.h"

namespace apps {
namespace hydro {

//! the 2d mesh type
using mesh_t = mesh__<2>;

// type aliases
using flux_data_t = flux_data__<mesh_t::num_dimensions>;

} // namespace 
} // namespace

