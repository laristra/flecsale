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
#include "flecsi/io/io_base.h"
#include "ale/mesh/burton/burton_mesh.h"
#include "ale/utils/errors.h"
#include "ale/utils/vtk.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
/// \class io_vtk_t io_vtk.h
/// \brief io_vtk_t provides a derived type of io_base.h and registrations
///   of the vtk file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_vtk_t : public flecsi::io_base_t<burton_mesh_t> {

  //! Default constructor
  burton_io_vtk_t() {}

  //============================================================================
  //! Implementation of vtk mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtk error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_t &m) 
  {

    std::cout << "Writing mesh to: " << name << std::endl;

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------
    // collect field data
    //--------------------------------------------------------------------------


    // collect the variables
    vector<string> variables;

    //----------------------------------------------------------------------------
    // nodal field data

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

    //--------------------------------------------------------------------------
    // WRITE HEADERS

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

    //--------------------------------------------------------------------------
    // WRITE DATA
    //--------------------------------------------------------------------------


    //----------------------------------------------------------------------------
    // coordinates

    // a temporary for vector values
    vector<vtk_real_t> vals( num_nodes * 3 );

    // get the coordinates from the mesh. unstructured always 3d
    for (auto v : m.vertices()) {
      const auto & coord = v->coordinates();
      for ( auto i=0; i<num_dims; i++ ) 
        vals[ 3*v.id() + i ] = coord[i];
      for ( auto i=num_dims; i<3; i++ )
        vals[ 3*v.id() + i ] = 0.0;
    } // for

    // write the coordinates to the file
    status = writer.write_points( vals, num_nodes, 3 );
    assert( status == 0 && "error with points" );

    //----------------------------------------------------------------------------
    // ellement connectivity

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
    assert( status == 0 && "error with element conn" );

    //----------------------------------------------------------------------------
    // nodal field data

    status = writer.start_point_data( num_nodes );
    assert( status == 0 && "error with point start" );

    vals.resize( num_nodes );
    vector<vtk_int_t> ivals( num_nodes );

    // node field buffer
    for(auto sf: rspav) {
      for(auto v: m.vertices()) vals[v.id()] = sf[v];
      status = writer.write_field( sf.label().c_str(), vals );
      assert( status == 0 && "error with point data" );
    } // for
    for(auto sf: ispav) {
      // cast int fields to real_t
      for(auto v: m.vertices()) ivals[v.id()] = sf[v];
      status = writer.write_field( sf.label().c_str(), ivals );
      assert( status == 0 && "error with point data" );
    } // for

    vals.resize( num_nodes * num_dims );
    for(auto vf: rvpav) {
      for(auto v: m.vertices()) {
        const auto & vec = vf[v];
        for ( auto i=0; i<num_dims; i++ ) 
          vals[ num_dims*v.id() + i ] = vec[i];
      } // for
      status = writer.write_field( vf.label().c_str(), vals, num_dims );
      assert( status == 0 && "error with cell data" );
    } // for

    //----------------------------------------------------------------------------
    // cell field data

    status = writer.start_cell_data( num_elem );
    assert( status == 0 && "error with cell start" );

    vals.resize( num_elem );
    ivals.resize( num_elem );

    // element field buffer
    for(auto sf: rspac) {
      for(auto c: m.cells()) vals[c.id()] = sf[c];
      status = writer.write_field( sf.label().c_str(), vals );
      assert( status == 0 && "error with cell data" );
    } // for
    for(auto sf: ispac) {
      // cast int fields to real_t
      for(auto c: m.cells()) ivals[c.id()] = sf[c];
      status = writer.write_field( sf.label().c_str(), ivals );
      assert( status == 0 && "error with cell data" );
    } // for
    
    // vectors
    vals.resize( num_elem * num_dims );
    for(auto vf: rvpac) {
      for(auto c: m.cells()) {
        const auto & vec = vf[c];
        for ( auto i=0; i<num_dims; i++ ) 
          vals[ num_dims*c.id() + i ] = vec[i];
      } // for
      status = writer.write_field( vf.label().c_str(), vals, num_dims );
      assert( status == 0 && "error with cell data" );
    } // for

    //--------------------------------------------------------------------------
    // CLOSE THE FILE
    //--------------------------------------------------------------------------

    // Create ZONE header
    status = writer.close();
    assert( status == 0 && "error with close" );


  } // io_vtk_t::write

  //============================================================================
  //! Implementation of vtk mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtk error code. 0 on success.
  //!
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_t &m) 
  {
    raise_implemented_error( "No vtk read functionality has been implemented" );
  };

}; // struct io_vtk_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_vtk_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_vtk_t.
//!
//! \return Pointer to io_base_t base class of io_vtk_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_t> * create_io_vtk()
{
  return new burton_io_vtk_t;
} // create_io_vtk


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtk" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_vtk_dat_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType(
    "vtk", create_io_vtk );


} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
