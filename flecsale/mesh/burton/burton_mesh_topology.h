/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines the burton mesh topology from the FleCSI topology type.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user incldues
#include "flecsale/mesh/burton/burton_types.h"
#include "flecsi/topology/mesh_types.h"


namespace flecsale {
namespace mesh {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
//! \brief Type for storing instance of template specialized low level mesh.
//! \tparam [in]  N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template < std::size_t N, bool Extra_Elements >
using burton_mesh_topology_t = 
  flecsi::topology::mesh_topology_t< burton_types_t<N, Extra_Elements> >;


} // namespace
} // namespace
} // namespace
