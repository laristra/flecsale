/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsi/io/io.h"
#include "mesh/burton/burton_mesh.h"



////////////////////////////////////////////////////////////////////////////////
// General Interface
////////////////////////////////////////////////////////////////////////////////

  
template< typename E >
auto filter_boundary( E && entities )  
{
  return 
    std::forward<decltype(entities)>(entities).filter( 
      [](auto e) { return e->is_boundary(); } 
    );
};


////////////////////////////////////////////////////////////////////////////////
// Alias mesh types
////////////////////////////////////////////////////////////////////////////////

namespace ale {
namespace mesh {

//! \brief The final 2d mesh type
using burton_mesh_2d_t = burton_mesh_t<2>;
//! \brief The final 3d mesh type
using burton_mesh_3d_t = burton_mesh_t<3>;

}
}

//! \brief Expose attributes and attachement sites to all namspaces.
//! This is horrible but it has to be done other wise users need to 
//! write stuff like flecsi_has_attribute_at( ale::mesh::burton::persistent,
//! ale::mesh::burton::vertex ).
using namespace ale::mesh::burton;

////////////////////////////////////////////////////////////////////////////////
// Delayed includes
////////////////////////////////////////////////////////////////////////////////

#include "mesh/burton/burton_io_exodus.h"
#include "mesh/burton/burton_io_tecplot.h"
#include "mesh/burton/burton_io_vtk.h"
#include "mesh/burton/burton_io_vtu.h"
#include "mesh/burton/burton_io_vtm.h"

////////////////////////////////////////////////////////////////////////////////
// load some things
////////////////////////////////////////////////////////////////////////////////

namespace ale {
namespace mesh {

//! \brief bring write/read mesh into the ale::mesh namespace
using flecsi::io::write_mesh;
//! \brief bring write/read mesh into the ale::mesh namespace
using flecsi::io::read_mesh;

}
}
