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
 * \file io_vtm.h
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
#  include <vtkCellType.h>
#  include <vtkMultiBlockDataSet.h>
#  include <vtkPoints.h>
#  include <vtkSmartPointer.h>
#  include <vtkUnstructuredGrid.h>
#  include <vtkXMLMultiBlockDataWriter.h>
#endif

#ifdef _DOUBLE_PRECISION_
#  undef _DOUBLE_PRECISION_
#  define DOUBLE_PRECISION
#endif

//! user includes
#include "flecsi/io/io_base.h"
#include "../../mesh/burton/burton_mesh.h"
#include "../../utils/errors.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
/// \class io_vtm_t io_vtm.h
/// \brief io_vtm_t provides a derived type of io_base.h and registrations
///   of the vtm file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_vtm_t : public flecsi::io_base_t<burton_mesh_t> {

  //! Default constructor
  burton_io_vtm_t() {}

  //============================================================================
  //! Implementation of vtm mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtm error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_t &m ) override
  {

#ifdef HAVE_VTK

    std::cout << "Writing mesh to: " << name << std::endl;

    // alias some types
    using std::string;
    using std::vector;

    using   mesh_t = burton_mesh_t;
    using   size_t = typename mesh_t::size_t;
    using   real_t = typename mesh_t::real_t;
    using integer_t= typename mesh_t::integer_t;
    using vector_t = typename mesh_t::vector_t;

    using vtkRealType = double;


    // get the regions
    auto region_cells = m.regions();
    auto num_regions = m.num_regions();

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    // vtkMultiBlockDataSet respresents multi-block datasets. See
    // the class documentation for more information.
    auto mb = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    
    for ( auto iblk=0; iblk<num_regions; iblk++ ) {

      // creat unstructured grid
      auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();


      // get points for this block
      std::map< mesh_t::vertex_t*, size_t > point_map;
      for ( auto c : region_cells[iblk] )
        for ( auto v : m.vertices(c) ) {
          auto id = point_map.size();
          point_map.emplace( std::make_pair(v,id) );
        }
      auto num_vertices = point_map.size();
      std::cout << num_vertices << std::endl;

      // create point storage
      auto points = vtkPoints::New();
      points->SetNumberOfPoints( num_vertices );
      
      // now create the vtk points, and create local ids
      for ( auto v : point_map ) {
        auto id = v.second;
        auto coord = v.first->coordinates();
        vtkRealType x[3];
        std::copy( coord.begin(), coord.end(), x );
        points->SetPoint( id, x );
      }

      // transfer the points
      ug->SetPoints(points);
      points->Delete();

      // create the cells
      for ( auto c : region_cells[iblk] ) {
        // get the vertices in this cell
        auto vs = m.vertices(c);
        auto n = vs.size();
        // copy them to the vtk type
        std::vector< vtkIdType > ids(n);
        std::transform( vs.begin(), vs.end(), ids.begin(),
                        [&](auto && v) { return point_map.at(v); } );
        // set the cell vertices
        ug->InsertNextCell(VTK_POLYGON, n, ids.data());
      }

      // Add the structured grid to the multi-block dataset
      mb->SetBlock(iblk, ug);
    
    } // block

    // write to file
    auto writer = vtkXMLMultiBlockDataWriter::New();
    writer->SetFileName( name.c_str() );
    writer->SetInputDataObject( mb );

    writer->Write();

    return 0;

#else

    std::cerr << "FLECSI not build with vtk support." << std::endl;
    std::exit(1);

    return -1;

#endif


  } // io_vtm_t::write

  //============================================================================
  //! Implementation of vtm mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtm error code. 0 on success.
  //!
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_t &m) override
  {
    raise_implemented_error( "No vtm read functionality has been implemented" );
  };

}; // struct io_vtm_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_vtm_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_vtm_t.
//!
//! \return Pointer to io_base_t base class of io_vtm_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_t> * create_io_vtm()
{
  return new burton_io_vtm_t;
} // create_io_vtm


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtm" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_vtm_dat_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType(
    "vtm", create_io_vtm );


} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
