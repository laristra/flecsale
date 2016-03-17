/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  // 
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
/*!
 * \file io_vtk.h
 * \date Initial file creation: Oct 07, 2015
 *
 ******************************************************************************/

#pragma once

//! system includes
#include <cstring>
#include <fstream>


//! user includes
#include "flecsi/io/io_vtk.h"
#include "ale/mesh/burton/burton_mesh.h"
#include "ale/utils/errors.h"
#include "ale/utils/vtk.h"


namespace flecsi {

//! bring burton mesh to front of namespace
using ale::mesh::burton_mesh_t;

//! Register file extension "plt" with factory.
static bool burton_vtk_dat_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType("vtk",
    flecsi::create_io_vtk<burton_mesh_t>);


////////////////////////////////////////////////////////////////////////////////
//! Implementation of vtk mesh write for burton specialization.
//!
//! \param[in] name Write burton mesh \e m to \e name.
//! \param[in] m Burton mesh to write to \e name.
//!
//! \return vtk error code. 0 on success.
//!
//! FIXME: should allow for const mesh_t &
////////////////////////////////////////////////////////////////////////////////
template<>
inline
int32_t 
flecsi::io_vtk_t<burton_mesh_t>::write( const std::string &name,
                                                   burton_mesh_t &m) 
{

  std::cout << "Writing mesh to: " << name << std::endl;

  //============================================================================
  // setup
  //============================================================================

  // alias some types
  using std::string;
  using std::vector;

  using   mesh_t = burton_mesh_t;
  using   real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t;
  using vtk_real_t = real_t;
  using vtk_int_t = int;

  using vtk_writer = ale::utils::vtk_writer;

  // get the general statistics
  auto num_dims  = m.num_dimensions();
  auto num_nodes = m.num_vertices();
  auto num_elem  = m.num_cells();
  constexpr auto num_nodes_per_elem = 4;

  // set the time
  auto soln_time = m.get_time();

  // variable extension for vectors
  std::string var_ext[3];
  var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

  // collect the variables
  vector<string> coordinates;

  // write the coordinate names
  coordinates.emplace_back( "x" );
  coordinates.emplace_back( "y" );

  //============================================================================
  // collect field data
  //============================================================================


  // collect the variables
  vector<string> variables;

  //----------------------------------------------------------------------------
  // nodal field data
  //----------------------------------------------------------------------------

  // real scalars persistent at vertices
  auto rspav = access_type_if(m, real_t, is_persistent_at(vertices));
  // int scalars persistent at vertices
  auto ispav = access_type_if(m, int, is_persistent_at(vertices));
  // real vectors persistent at vertices
  auto rvpav = access_type_if(m, vector_t, is_persistent_at(vertices));

  // fill node variable names array
  for(auto sf: rspav) 
    variables.emplace_back( sf.label() );
  for(auto sf: ispav) 
    variables.emplace_back( sf.label() );
  for(auto vf: rvpav) {
    auto label = vf.label();
    for(int d=0; d < num_dims; ++d) {
      auto dim_label = label + var_ext[d];
      variables.emplace_back( dim_label );
    } // for
  } // for
  

  //----------------------------------------------------------------------------
  // element field data
  //----------------------------------------------------------------------------

  // real scalars persistent at cells
  auto rspac = access_type_if(m, real_t, is_persistent_at(cells));
  // int scalars persistent at cells
  auto ispac = access_type_if(m, int, is_persistent_at(cells));
  // real vectors persistent at cells
  auto rvpac = access_type_if(m, vector_t, is_persistent_at(cells));


  // fill element variable names array
  for(auto sf: rspac) 
    variables.emplace_back( sf.label() );
  for(auto sf: ispac) 
    variables.emplace_back( sf.label() );
  for(auto vf: rvpac) {
    auto label = vf.label();
    for(int d=0; d < num_dims; ++d) {
      auto dim_label = label + var_ext[d];
      variables.emplace_back( dim_label );
    } // for
  } // for

  //============================================================================
  // WRITE HEADERS
  //============================================================================

  // create the variable name string
  string var_string;
  for ( const auto & label : variables )
    var_string += label + " ";

  // open the file
  vtk_writer writer;
  auto status = writer.open( name.c_str() );  
  assert( status == 0 && "error with open" );
  
  // Create ZONE header
  status = writer.init( "Quadrilateral Zone" );
  assert( status == 0 && "error with header" );

  //============================================================================
  // WRITE DATA
  //============================================================================


  //----------------------------------------------------------------------------
  // coordinates
  //----------------------------------------------------------------------------

  // a temporary for vector values
  vector<vtk_real_t> vals( num_nodes * num_dims );

  // get the coordinates from the mesh.
  for (auto v : m.vertices()) {
    vals[ 0*num_nodes + v.id() ] = v->coordinates()[0];
    vals[ 1*num_nodes + v.id() ] = v->coordinates()[1];
  } // for

  // write the coordinates to the file
  status = writer.write_points( vals, num_nodes, num_dims );


  //----------------------------------------------------------------------------
  // ellement connectivity
  //----------------------------------------------------------------------------

  // element definitions
  vector< vector<vtk_int_t> > elem_conn( num_elem );
  for (auto c : m.cells()) {
    auto cid = c.id();
    elem_conn[cid].resize( num_nodes_per_elem );
    auto vid = 0;
    for (auto v : m.vertices(c)) {
      elem_conn[cid][vid++] = v.id();
    } // for
  } // for

  status = writer.write_elements( elem_conn, vtk_writer::cell_type_t::vtk_quad );
  assert( status == 0 && "error with TECNOD" );


  //============================================================================
  // CLOSE THE FILE
  //============================================================================

  // Create ZONE header
  status = writer.close();
  assert( status == 0 && "error with close" );


} // io_vtk_t::write

////////////////////////////////////////////////////////////////////////////////
//! Implementation of vtk mesh read for burton specialization.
//!
//! \param[in] name Read burton mesh \e m to \e name.
//! \param[in] m Burton mesh to Read to \e name.
//!
//! \return vtk error code. 0 on success.
//!
////////////////////////////////////////////////////////////////////////////////
template<>
inline
int32_t 
flecsi::io_vtk_t<burton_mesh_t>::read( const std::string &name,
                                                  burton_mesh_t &m) 
{
  raise_implemented_error( "No vtk read functionality has been implemented" );
};

} // namespace flecsi

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
