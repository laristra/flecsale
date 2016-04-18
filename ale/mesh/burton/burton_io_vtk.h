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

// vtk doesnt like double-precision
#ifdef DOUBLE_PRECISION
#  undef DOUBLE_PRECISION
#  define _DOUBLE_PRECISION_
#endif

#ifdef HAVE_VTK
#  include <vtkUnstructuredGridReader.h>
#  include <vtkUnstructuredGridWriter.h>
#endif


#ifdef _DOUBLE_PRECISION_
#  undef _DOUBLE_PRECISION_
#  define DOUBLE_PRECISION
#endif


//! user includes
#include "flecsi/io/io_base.h"
#include "ale/mesh/burton/burton_mesh.h"
#include "ale/mesh/vtk_utils.h"
#include "ale/utils/errors.h"


namespace ale {
namespace mesh {


////////////////////////////////////////////////////////////////////////////////
/// \class io_vtk_t io_vtk.h
/// \brief io_vtk_t provides a derived type of io_base.h and registrations
///   of the vtk file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_vtk_t : public flecsi::io_base_t<burton_mesh_2d_t> {

  //! Default constructor
  burton_io_vtk_t() {}

#if HAVE_VTK


  //============================================================================
  //! Implementation of vtu mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtu error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //!
  //! \remark this uses in vtk library writer
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_2d_t &m, bool binary ) override
  {

    std::cout << "Writing mesh to: " << name << std::endl;

    // convert to vtk
    auto ug = to_vtk( m );

    // write to file
    auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName( name.c_str() );
    writer->SetInputDataObject( ug );
    if ( binary ) writer->SetFileTypeToBinary();
    else          writer->SetFileTypeToASCII();

    writer->Write();

    return 0;


  } // io_vtu_t::write


  //============================================================================
  //! Implementation of vtu mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtu error code. 0 on success.
  //!
  //! \remark this uses in vtk library reader
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_2d_t &m) override
  {

    std::cout << "Reading mesh from: " << name << std::endl;


    // Read solution
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName( name.c_str() );
    reader->Update();
    auto ug = reader->GetOutput();
    
    // convert vtk solution to a mesh
    m = to_mesh<burton_mesh_2d_t>( ug );

    return 0;

  };


#else

  //============================================================================
  //! Implementation of vtk mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtk error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //! 
  //! \remark this uses in house vtk writer
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_2d_t &m, bool binary ) override
  {

    std::cout << "Writing mesh to: " << name << std::endl;

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    // alias some types
    using std::vector;

    using   mesh_t = burton_mesh_2d_t;
    using   size_t = typename mesh_t::size_t;
    using   real_t = typename mesh_t::real_t;
    using integer_t= typename mesh_t::integer_t;
    using vector_t = typename mesh_t::vector_t;


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    // get the general statistics
    auto num_dims  = m.num_dimensions();
    auto num_nodes = m.num_vertices();
    auto num_elem  = m.num_cells();

    // set the time
    auto soln_time = m.time();

    //--------------------------------------------------------------------------
    // WRITE HEADERS
    //--------------------------------------------------------------------------
    
    using ale::utils::vtk_writer;

    // open the file
    vtk_writer writer;
    auto status = writer.open( name.c_str(), binary );  
    assert( status == 0 && "error with open" );
  
    // Create ZONE header
    status = writer.init( "Mesh Zone" );
    assert( status == 0 && "error with header" );

    //--------------------------------------------------------------------------
    // WRITE DATA
    //--------------------------------------------------------------------------


    //----------------------------------------------------------------------------
    // coordinates

    // a temporary for vector values
    vector<real_t> vals( num_nodes * 3 );

    // get the coordinates from the mesh. unstructured always 3d
    for (auto v : m.vertices()) {
      const auto & coord = v->coordinates();
      auto base = 3*v.id();
      for ( auto i=0; i<num_dims; i++ ) vals[ base + i ] = coord[i];
      for ( auto i=num_dims; i<3; i++ ) vals[ base + i ] = 0.0;
    } // for

    // write the coordinates to the file
    status = writer.write_points( vals, num_nodes, 3 );
    assert( status == 0 && "error with points" );

    //----------------------------------------------------------------------------
    // ellement connectivity
    
    vector< vector<integer_t> > elem_conn( num_elem );
    vector< vtk_writer::cell_type_t > elem_types( num_elem );

    // element definitions
    for (auto c : m.cells()) {
      auto cid = c.id();
      auto elem_verts = m.vertices(c);
      auto num_nodes_per_elem = elem_verts.size();
      elem_conn[cid].resize( num_nodes_per_elem );
      size_t vid = 0;
      for (auto v : elem_verts) elem_conn[cid][vid++] = v.id();
      elem_types[cid] = vtk_writer::cell_type_t::polygon;
    } // for

    status = writer.write_elements( elem_conn, elem_types.data() );
    assert( status == 0 && "error with element conn" );

    //----------------------------------------------------------------------------
    // nodal field data

    status = writer.start_point_data( num_nodes );
    assert( status == 0 && "error with point start" );

    vals.resize( num_nodes );
    vector<integer_t> ivals( num_nodes );

    // real scalars persistent at vertices
    auto rspav = access_type_if(m, real_t, is_persistent_at(m,vertices));
    for(auto sf: rspav) {
      auto label = validate_string( sf.label() );
      for(auto v: m.vertices()) vals[v.id()] = sf[v];
      status = writer.write_field( label.c_str(), vals );
      assert( status == 0 && "error with point data" );
    } // for

    // int scalars persistent at vertices
    auto ispav = access_type_if(m, integer_t, is_persistent_at(m,vertices));
    for(auto sf: ispav) {
      auto label = validate_string( sf.label() );
      for(auto v: m.vertices()) ivals[v.id()] = sf[v];
      status = writer.write_field( label.c_str(), ivals );
      assert( status == 0 && "error with point data" );
    } // for

    vals.resize( num_nodes * num_dims );

    // real vectors persistent at vertices
    auto rvpav = access_type_if(m, vector_t, is_persistent_at(m,vertices));
    for(auto vf: rvpav) {
      auto label = validate_string( vf.label() );
      for(auto v: m.vertices()) {
        const auto & vec = vf[v];
        for ( auto i=0; i<num_dims; i++ ) 
          vals[ num_dims*v.id() + i ] = vec[i];
      } // for
      status = writer.write_field( label.c_str(), vals, num_dims );
      assert( status == 0 && "error with cell data" );
    } // for

    //----------------------------------------------------------------------------
    // cell field data

    status = writer.start_cell_data( num_elem );
    assert( status == 0 && "error with cell start" );

    vals.resize( num_elem );
    ivals.resize( num_elem );

    //----------------------------------------------------------------------------
    // element field data

    // real scalars persistent at cells
    auto rspac = access_type_if(m, real_t, is_persistent_at(m,cells));
    for(auto sf: rspac) {
      auto label = validate_string( sf.label() );
      for(auto c: m.cells()) vals[c.id()] = sf[c];
      status = writer.write_field( label.c_str(), vals );
      assert( status == 0 && "error with cell data" );
    } // for

    // int scalars persistent at cells
    auto ispac = access_type_if(m, integer_t, is_persistent_at(m,cells));
    for(auto sf: ispac) {
      auto label = validate_string( sf.label() );
      for(auto c: m.cells()) ivals[c.id()] = sf[c];
      status = writer.write_field( label.c_str(), ivals );
      assert( status == 0 && "error with cell data" );
    } // for
    
    vals.resize( num_elem * num_dims );

    // real vectors persistent at cells
    auto rvpac = access_type_if(m, vector_t, is_persistent_at(m,cells));
    for(auto vf: rvpac) {
      auto label = validate_string( vf.label() );
      for(auto c: m.cells()) {
        const auto & vec = vf[c];
        for ( auto i=0; i<num_dims; i++ ) 
          vals[ num_dims*c.id() + i ] = vec[i];
      } // for
      status = writer.write_field( label.c_str(), vals, num_dims );
      assert( status == 0 && "error with cell data" );
    } // for

    //--------------------------------------------------------------------------
    // CLOSE THE FILE
    //--------------------------------------------------------------------------

    // Create ZONE header
    status = writer.close();
    assert( status == 0 && "error with close" );

    return status;

  } // io_vtk_t::write


  //============================================================================
  //! Implementation of vtk mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtk error code. 0 on success.
  //!
  //! \remark no read support using in house vtk writer
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_2d_t &m) override
  {
    raise_implemented_error( "No vtk read functionality has been implemented" );
  };

#endif // HAVE_VTK

  //============================================================================
  //! Implementation of vtk mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtk error code. 0 on success.
  //!
  //! \remark default to ascii
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_2d_t &m ) override
  {
    write( name, m, true );
  }

}; // struct io_vtk_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_vtk_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_vtk_t.
//!
//! \return Pointer to io_base_t base class of io_vtk_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_2d_t> * create_io_vtk()
{
  return new burton_io_vtk_t;
} // create_io_vtk


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtk" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_vtk_dat_registered =
  flecsi::io_factory_t<burton_mesh_2d_t>::instance().registerType(
    "vtk", create_io_vtk );


} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
