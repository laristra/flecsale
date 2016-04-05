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

#ifdef HAVE_VTK
#  include <vtkCellArray.h>
#  include <vtkCellData.h>
#  include <vtkCellType.h>
#  include <vtkPointData.h>
#  include <vtkPoints.h>
#  include <vtkSmartPointer.h>
#  include <vtkUnstructuredGrid.h>
#endif

//! user includes
#include "../../utils/vtk.h"


namespace ale {
namespace mesh {


////////////////////////////////////////////////////////////////////////////////
//! \brief convert a burton mesh to a vtk unstructured grid
////////////////////////////////////////////////////////////////////////////////
static auto to_vtk( burton_mesh_t & m ) 
{

  // alias some types
  using std::vector;

  using   mesh_t = burton_mesh_t;
  using   real_t = typename mesh_t::real_t;
  using integer_t= typename mesh_t::integer_t;
  using vector_t = typename mesh_t::vector_t;

  using utils::vtk_array_t;

  // a lambda function for validating strings
  auto validate_string = []( auto && str ) {
    return std::forward<decltype(str)>(str);
  };

  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------

  // some general mesh stats
  auto num_dims = m.num_dimensions();
  auto num_vertices = m.num_vertices();
  auto num_cells = m.num_cells();
    
  // creat unstructured grid
  auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // create points
  auto points = vtkPoints::New();
  points->SetNumberOfPoints( num_vertices );
  for ( auto v : m.vertices() ) {
    auto id = v.id();
    auto & coord = v->coordinates();
    real_t x[3] = {0, 0, 0};
    std::copy( coord.begin(), coord.end(), x );
    points->SetPoint( id, x );
  }

  // transfer the points
  ug->SetPoints(points);
  points->Delete();

  // create the cells
  for ( auto c : m.cells() ) {
    // get the vertices in this cell
    auto vs = m.vertices(c);
    auto n = vs.size();
    // copy them to the vtk type
    vector< vtkIdType > ids(n);
    std::transform( vs.begin(), vs.end(), ids.begin(),
                    [](auto && v) { return v.id(); } );
    // set the cell vertices
    ug->InsertNextCell(VTK_POLYGON, n, ids.data());
  }
    

  //----------------------------------------------------------------------------
  // write field data
  //----------------------------------------------------------------------------
    
  //----------------------------------------------------------------------------
  // nodal field data header

  // get the point data object
  auto pd = ug->GetPointData();

  // real scalars persistent at vertices
  auto rvals = vtk_array_t<real_t>::type::New();
  rvals->SetNumberOfValues( num_vertices );

  auto rspav = access_type_if(m, real_t, is_persistent_at(vertices));
  for(auto sf: rspav) {
    auto label = validate_string( sf.label() );      
    rvals->SetName( label.c_str() );
    for(auto v: m.vertices()) rvals->SetValue( v.id(), sf[v] );
    pd->AddArray( rvals );
  } // for

  rvals->Delete();

  // int scalars persistent at vertices
  auto ivals = vtk_array_t<integer_t>::type::New();
  ivals->SetNumberOfValues( num_vertices );

  auto ispav = access_type_if(m, integer_t, is_persistent_at(vertices));
  for(auto sf: ispav) {
    auto label = validate_string( sf.label() );
    ivals->SetName( label.c_str() );
    for(auto v: m.vertices()) ivals->SetValue( v.id(), sf[v] );
    pd->AddArray( ivals );
  } // for
    
  ivals->Delete();

  // real vectors persistent at vertices
  auto vvals = vtk_array_t<real_t>::type::New();
  vvals->SetNumberOfComponents( 3 ); // always 3d
  vvals->SetNumberOfTuples( num_vertices );

  auto rvpav = access_type_if(m, vector_t, is_persistent_at(vertices));
  for(auto vf: rvpav) {
    auto label = validate_string( vf.label() );
    vvals->SetName( label.c_str() );
    for(auto v: m.vertices()) {
      real_t to_vals[3] = {0, 0, 0};
      auto & from_vals = vf[v];
      std::copy( from_vals.begin(), from_vals.end(), to_vals );
      vvals->SetTuple( v.id(), to_vals );
    } // for
    pd->AddArray( vvals );
  } // for

  vvals->Delete();

  //----------------------------------------------------------------------------
  // cell field data header

  // get the point data object
  auto cd = ug->GetCellData();


  // real scalars persistent at cells
  rvals = vtk_array_t<real_t>::type::New();
  rvals->SetNumberOfValues( num_cells );

  auto rspac = access_type_if(m, real_t, is_persistent_at(cells));
  for(auto sf: rspac) {
    auto label = validate_string( sf.label() );      
    rvals->SetName( label.c_str() );
    for(auto c: m.cells()) rvals->SetValue( c.id(), sf[c] );
    cd->AddArray( rvals );
  } // for

  rvals->Delete();

  // int scalars persistent at cells
  ivals = vtk_array_t<integer_t>::type::New();
  ivals->SetNumberOfValues( num_cells );

  auto ispac = access_type_if(m, integer_t, is_persistent_at(cells));
  for(auto sf: ispac) {
    auto label = validate_string( sf.label() );
    ivals->SetName( label.c_str() );
    for(auto c: m.cells()) ivals->SetValue( c.id(), sf[c] );
    cd->AddArray( ivals );
  } // for
    
  ivals->Delete();

  // real vectors persistent at cells
  vvals = vtk_array_t<real_t>::type::New();
  vvals->SetNumberOfComponents( 3 ); // always 3d
  vvals->SetNumberOfTuples( num_cells );

  auto rvpac = access_type_if(m, vector_t, is_persistent_at(cells));
  for(auto vf: rvpac) {
    auto label = validate_string( vf.label() );
    vvals->SetName( label.c_str() );
    for(auto c: m.cells()) {
      real_t to_vals[3] = {0, 0, 0};
      auto & from_vals = vf[c];
      std::copy( from_vals.begin(), from_vals.end(), to_vals );
      vvals->SetTuple( c.id(), to_vals );
    } // for
    cd->AddArray( vvals );
  } // for

  vvals->Delete();

  //----------------------------------------------------------------------------
  // return mesh
  //----------------------------------------------------------------------------
  return ug;

}

////////////////////////////////////////////////////////////////////////////////
//! \brief convert a vtk unstructured grid to a burton mesh
////////////////////////////////////////////////////////////////////////////////
static void to_mesh( vtkUnstructuredGrid* ug, burton_mesh_t & m ) 
{


  // alias some types
  using std::vector;

  using   mesh_t = burton_mesh_t;
  using   real_t = typename mesh_t::real_t;
  using integer_t= typename mesh_t::integer_t;
  using vector_t = typename mesh_t::vector_t;
  using  point_t = typename mesh_t::point_t;
  using vertex_t = typename mesh_t::vertex_t;

  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------
  
  // some general mesh stats
  auto num_dims = m.num_dimensions();


  //----------------------------------------------------------------------------
  // Points
  //----------------------------------------------------------------------------
    
  // get points
  auto points = ug->GetPoints();
  auto num_vertices = points->GetNumberOfPoints();

  m.init_parameters(num_vertices);

  // create points
  vector<vertex_t *> vs;
  vs.reserve( num_vertices );

  for (size_t i = 0; i < num_vertices; ++i) {
    real_t x[3] = {0, 0, 0};
    points->GetPoint( i, x );
    point_t p = { static_cast<real_t>(x[0]), 
                  static_cast<real_t>(x[1]) };
    auto v = m.create_vertex(p);
    vs.emplace_back( std::move(v) );
  } // for

  //----------------------------------------------------------------------------
  // Cells
  //----------------------------------------------------------------------------
  
  // get cells
  auto cells = ug->GetCells();
  auto num_cells = cells->GetNumberOfCells();

  vector<vertex_t *> elem_vs;
  vtkIdType npts, *pts;

  // create the cells
  cells->InitTraversal();
  while( cells->GetNextCell(npts, pts) ) {
    // reset storage
    elem_vs.clear();
    elem_vs.reserve( npts );
    // get the points
    for ( auto v=0;  v<npts; v++ )
      elem_vs.emplace_back( vs[ pts[v] ] );
    // create acual cell
    auto c = m.create_cell( elem_vs );                
  }
    

  //----------------------------------------------------------------------------
  // Finish up
  //----------------------------------------------------------------------------

  m.init();

}

} // namespace
} // namespace
