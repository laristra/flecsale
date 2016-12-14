/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Defines a test fixture.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// test include
#include "burton_test_base.h"

////////////////////////////////////////////////////////////////////////////////
//! \brief A utility for creating data
////////////////////////////////////////////////////////////////////////////////
template< typename M >
void create_data(M & mesh ) 
{

  //! \brief the mesh float type
  using real_t   = typename M::real_t;
  //! \brief the mesh int type
  using integer_t= typename M::integer_t;
  //! \brief the vector type
  using vector_t = typename M::vector_t;

  // register
  register_data(mesh, hydro, pressure, real_t, dense, 1, cells);
  register_data(mesh, hydro, region, integer_t, dense, 1, cells);
  register_data(mesh, hydro, velocity, vector_t, dense, 1, vertices);
  // access
  auto p = get_accessor(mesh, hydro, pressure, real_t, dense, 0);
  auto r = get_accessor(mesh, hydro, region, integer_t, dense, 0);
  auto velocity = get_accessor(mesh, hydro, velocity, vector_t, dense, 0);
  // set attributes
  p.attributes().set(persistent);
  r.attributes().set(persistent);
  velocity.attributes().set(persistent);
  // initialize
  for(auto c: mesh.cells()) {
    p[c] = c.id();
    r[c] = mesh.num_cells() - c.id();
  } // for
    // vertices
  for (auto v: mesh.vertices()) {
    velocity[v][0] = v.id();
    velocity[v][1] = 2.0*v.id();
  } // for
}



////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh
////////////////////////////////////////////////////////////////////////////////
class burton_io : public burton_test_base {};
