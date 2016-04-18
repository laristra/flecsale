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
 * \file io_vtu.h
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
#  include <vtkXMLUnstructuredGridReader.h>
#  include <vtkXMLUnstructuredGridWriter.h>
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
/// \class io_vtu_t io_vtu.h
/// \brief io_vtu_t provides a derived type of io_base.h and registrations
///   of the vtu file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_vtu_t : public flecsi::io_base_t<burton_mesh_2d_t> {

  //! Default constructor
  burton_io_vtu_t() {}

  //============================================================================
  //! Implementation of vtu mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtu error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_2d_t &m ) override
  {

#ifdef HAVE_VTK

    std::cout << "Writing mesh to: " << name << std::endl;

    // convert to vtk
    auto ug = to_vtk( m );

    // write to file
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName( name.c_str() );
    writer->SetInputDataObject( ug );

    writer->Write();

    return 0;

#else

    raise_implemented_error( "FLECSI not build with vtk support." );

    return -1;

#endif


  } // io_vtu_t::write


  //============================================================================
  //! Implementation of vtu mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtu error code. 0 on success.
  //!
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_2d_t &m) override
  {
#ifdef HAVE_VTK

    std::cout << "Reading mesh from: " << name << std::endl;

    // Read solution
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName( name.c_str() );
    reader->Update();
    auto ug = reader->GetOutput();

    // convert vtk solution to a mesh
    m = to_mesh<burton_mesh_2d_t>( ug );


    return 0;

#else

    raise_implemented_error( "FLECSI not build with vtk support." );

    return -1;

#endif

  };

}; // struct io_vtu_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_vtu_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_vtu_t.
//!
//! \return Pointer to io_base_t base class of io_vtu_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_2d_t> * create_io_vtu()
{
  return new burton_io_vtu_t;
} // create_io_vtu


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtu" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_vtu_dat_registered =
  flecsi::io_factory_t<burton_mesh_2d_t>::instance().registerType(
    "vtu", create_io_vtu );


} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
