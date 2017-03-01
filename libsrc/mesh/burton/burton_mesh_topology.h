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
#include "mesh/burton/burton_types.h"
#include "flecsi/topology/mesh_types.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \brief Type for storing instance of template specialized low level mesh.
//! \tparam [in]  N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template < std::size_t N >
using burton_mesh_topology_t = 
  flecsi::topology::mesh_topology_t< burton_types_t<N> >;


//! Two dimensional specialization of the mesh topology.
using burton_2d_mesh_topology_t = burton_mesh_topology_t<2>; 
//! Three dimensional specialization of the mesh topology.
using burton_3d_mesh_topology_t = burton_mesh_topology_t<3>; 


} // namespace
} // namespace
