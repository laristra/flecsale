/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This file defines the vtk reader and writer.  
/// \details There are two versions, one uses the vtk library and the other
///          uses a native reader.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsi/io/io_base.h"
#include "mesh/burton/burton_mesh.h"
#include "mesh/vtk_utils.h"
#include "utils/errors.h"


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

// system includes
#include <cstring>
#include <fstream>


namespace ale {
namespace mesh {


////////////////////////////////////////////////////////////////////////////////
/// \brief This is the mesh reader and writer based on the vtk format.
/// \tparam N The number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_io_vtk_t : public flecsi::io::io_base_t< burton_mesh_t<N> > {

public:

  //! the mesh type
  using mesh_t = burton_mesh_t<N>;

  //! Default constructor
  burton_io_vtk_t() {}

#if HAVE_VTK


  //============================================================================
  //! \brief Implementation of vtk mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //!
  //! \remark this uses in vtk library writer
  //============================================================================
  int write( const std::string &name, mesh_t &m, bool binary ) override
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


  }


  //============================================================================
  //! \brief Implementation of vtk mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return error code. 0 on success.
  //!
  //! \remark this uses in vtk library reader
  //============================================================================
  int read( const std::string &name, mesh_t &m) override
  {

    std::cout << "Reading mesh from: " << name << std::endl;


    // Read solution
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName( name.c_str() );
    reader->Update();
    auto ug = reader->GetOutput();
    
    // convert vtk solution to a mesh
    m = to_mesh<mesh_t>( ug );

    return 0;

  };


#else

  //============================================================================
  //! \brief Implementation of vtk mesh write for burton specialization.
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
  int write( const std::string &name, mesh_t &m, bool binary ) override
  {

    std::cout << "Writing mesh to: " << name << std::endl;

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    // alias some types
    using std::vector;

    using   size_t = typename mesh_t::size_t;
    using   real_t = typename mesh_t::real_t;
    using integer_t= typename mesh_t::integer_t;
    using vector_t = typename mesh_t::vector_t;


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    // get the general statistics
    constexpr auto num_dims  = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();
    auto num_elem  = m.num_cells();

    // set the time
    auto soln_time = m.time();

    //--------------------------------------------------------------------------
    // WRITE HEADERS
    //--------------------------------------------------------------------------
    
    using ale::io::vtk_writer;

    // open the file
    vtk_writer writer;
    auto status = writer.open( name.c_str(), binary );  
    assert( status == 0 && "error with open" );
  
    // Create ZONE header
    status = writer.init( "vtk output" );
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
      for ( int i=0; i<num_dims; i++ ) vals[ base + i ] = coord[i];
      for ( int i=num_dims; i<3; i++ ) vals[ base + i ] = 0.0;
    } // for

    // write the coordinates to the file
    status = writer.write_points( vals, num_nodes, 3 );
    assert( status == 0 && "error with points" );

    //----------------------------------------------------------------------------
    // ellement connectivity
    
    vector< vector<vtkIdType> > elem_conn( num_elem );
    vector< vtk_writer::cell_type_t > elem_types( num_elem );

    // element definitions
    if ( N ==2 ) {

      for (auto c : m.cells()) {
        auto cid = c.id();
        auto elem_verts = m.vertices(c);
        auto num_nodes_per_elem = elem_verts.size();
        elem_conn[cid].resize( num_nodes_per_elem );
        size_t vid = 0;
        for (auto v : elem_verts) elem_conn[cid][vid++] = v.id();
        elem_types[cid] = vtk_writer::cell_type_t::polygon;
      } // for

    }
    else if ( N == 3 ) {

      for (auto c : m.cells()) {
        auto cid = c.id();
        // get the element verts and faces
        auto elem_verts = m.vertices(c);
        auto elem_faces = m.faces(c);
        auto num_elem_verts = elem_verts.size();
        auto num_elem_faces = elem_faces.size();
        // get the total number of vertices
        auto tot_verts = std::accumulate( 
          elem_faces.begin(), elem_faces.end(), static_cast<size_t>(0),
          [&m](auto sum, auto f) { return sum + m.vertices(f).size(); }
        );
        // the list of faces that vtk requires contains the total number of faces
        // AND the number of points in each face AND the point ids themselves.
        elem_conn[cid].reserve( tot_verts + num_elem_faces + 1);
        elem_conn[cid].emplace_back( num_elem_faces );
        // now loop over each face and add the face data
        for ( auto f : elem_faces ) {
          auto face_elems = m.cells(f);
          auto face_verts = m.vertices(f);
          auto num_face_verts = face_verts.size();
          // copy the face vert ids
          vector< size_t > face_vert_ids( num_face_verts );
          std::transform( 
            face_verts.begin(), face_verts.end(), face_vert_ids.begin(),
            [](auto && v) { return v.id(); } 
          );
          // check the direction of the vertices
          if ( face_elems[0] != c ) 
            std::reverse( face_vert_ids.begin(), face_vert_ids.end() );
          // now copy them to the global array
          elem_conn[cid].emplace_back( num_face_verts );
          for ( auto v : face_vert_ids )
            elem_conn[cid].emplace_back( v );
        }
        elem_types[cid] = vtk_writer::cell_type_t::polyhedron;
      } // for

    } 
    else  {

      raise_logic_error( "No support for N="<<N<<" dimensions" );

    }

    status = writer.write_elements( elem_conn, elem_types.data() );
    assert( status == 0 && "error with element conn" );

    //----------------------------------------------------------------------------
    // nodal field data

    status = writer.start_point_data( num_nodes );
    assert( status == 0 && "error with point start" );

    vals.resize( num_nodes );
    vector<integer_t> ivals( num_nodes );

    // real scalars persistent at vertices
    auto rspav = flecsi_get_accessors_all(
      m, real_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
    );
    for(auto sf: rspav) {
      auto label = validate_string( sf.label() );
      for(auto v: m.vertices()) vals[v.id()] = sf[v];
      status = writer.write_field( label.c_str(), vals );
      assert( status == 0 && "error with point data" );
    } // for

    // int scalars persistent at vertices
    auto ispav = flecsi_get_accessors_all(
      m, integer_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
    );
    for(auto sf: ispav) {
      auto label = validate_string( sf.label() );
      for(auto v: m.vertices()) ivals[v.id()] = sf[v];
      status = writer.write_field( label.c_str(), ivals );
      assert( status == 0 && "error with point data" );
    } // for

    vals.resize( num_nodes * num_dims );

    // real vectors persistent at vertices
    auto rvpav = flecsi_get_accessors_all(
      m, vector_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
    );
    for(auto vf: rvpav) {
      auto label = validate_string( vf.label() );
      for(auto v: m.vertices()) {
        const auto & vec = vf[v];
        for ( int i=0; i<num_dims; i++ ) 
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
    auto rspac = flecsi_get_accessors_all(
      m, real_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
    );
    for(auto sf: rspac) {
      auto label = validate_string( sf.label() );
      for(auto c: m.cells()) vals[c.id()] = sf[c];
      status = writer.write_field( label.c_str(), vals );
      assert( status == 0 && "error with cell data" );
    } // for

    // int scalars persistent at cells
    auto ispac = flecsi_get_accessors_all(
      m, integer_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
    );
    for(auto sf: ispac) {
      auto label = validate_string( sf.label() );
      for(auto c: m.cells()) ivals[c.id()] = sf[c];
      status = writer.write_field( label.c_str(), ivals );
      assert( status == 0 && "error with cell data" );
    } // for
    
    vals.resize( num_elem * num_dims );

    // real vectors persistent at cells
    auto rvpac = flecsi_get_accessors_all(
      m, vector_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
    );
    for(auto vf: rvpac) {
      auto label = validate_string( vf.label() );
      for(auto c: m.cells()) {
        const auto & vec = vf[c];
        for ( int i=0; i<num_dims; i++ ) 
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
  //! \brief Implementation of vtk mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtk error code. 0 on success.
  //!
  //! \remark no read support using in house vtk writer
  //============================================================================
  int read( const std::string &name, mesh_t &m) override
  {
    raise_implemented_error( "No vtk read functionality has been implemented" );
  };

#endif // HAVE_VTK

  //============================================================================
  //! \brief Implementation of vtk mesh write for burton specialization.
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
  int write( const std::string &name, mesh_t &m ) override
  {
    return write( name, m, true );
  }

}; // struct io_vtk_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_vtk_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_vtk_t.
//!
//! \return Pointer to io_base_t base class of io_vtk_t.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
inline flecsi::io::io_base_t< burton_mesh_t<N> > * create_io_vtk()
{
  return new burton_io_vtk_t<N>;
} // create_io_vtk


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtk" with factory.
////////////////////////////////////////////////////////////////////////////////
//! @{
static bool burton_2d_vtk_registered =
  flecsi::io::io_factory_t< burton_mesh_t<2> >::instance().registerType(
    "vtk", create_io_vtk );

static bool burton_3d_vtk_registered =
  flecsi::io::io_factory_t< burton_mesh_t<3> >::instance().registerType(
    "vtk", create_io_vtk );
//! @}

} // namespace mesh
} // namespace ale
