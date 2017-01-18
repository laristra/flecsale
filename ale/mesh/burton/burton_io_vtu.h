/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Provides the functionality for writing meshes in vtu format.
////////////////////////////////////////////////////////////////////////////////

#pragma once


// user includes
#include "flecsi/io/io_base.h"
#include "ale/mesh/burton/burton_mesh.h"
#ifdef HAVE_VTK
#include "ale/mesh/vtk_utils.h"
#endif
#include "ale/utils/errors.h"


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

// system includes
#include <cstring>
#include <fstream>

namespace ale {
namespace mesh {



////////////////////////////////////////////////////////////////////////////////
/// \brief This is the mesh reader and writer based on the vtk format.
////////////////////////////////////////////////////////////////////////////////
template<std::size_t N>
class burton_io_vtu_t : public flecsi::io::io_base_t<burton_mesh_t<N>> {

public:

  //! the mesh type
  using mesh_t = burton_mesh_t<N>;

  //! Default constructor
  burton_io_vtu_t() {}

  //============================================================================
  //! \brief Implementation of vtu mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtu error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, mesh_t &m ) override
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
  int32_t read( const std::string &name, mesh_t &m) override
  {
#ifdef HAVE_VTK

    std::cout << "Reading mesh from: " << name << std::endl;

    // Read solution
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName( name.c_str() );
    reader->Update();
    auto ug = reader->GetOutput();

    // convert vtk solution to a mesh
    m = to_mesh<mesh_t>( ug );


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
template< std::size_t N >
inline flecsi::io::io_base_t<burton_mesh_t<N>> * create_io_vtu()
{
  return new burton_io_vtu_t<N>;
} // create_io_vtu


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtu" with factory.
////////////////////////////////////////////////////////////////////////////////
//! @{
static bool burton_2d_vtu_registered =
  flecsi::io::io_factory_t<burton_mesh_t<2>>::instance().registerType(
    "vtu", create_io_vtu );

static bool burton_3d_vtu_registered =
  flecsi::io::io_factory_t<burton_mesh_t<3>>::instance().registerType(
    "vtu", create_io_vtu );
//! @}

} // namespace mesh
} // namespace ale
