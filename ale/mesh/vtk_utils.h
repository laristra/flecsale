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

#  include <vtkDoubleArray.h>
#  include <vtkFloatArray.h>
#  include <vtkIntArray.h>
#  include <vtkLongArray.h>
#  include <vtkLongLongArray.h>
#endif

//! user includes
#include "ale/io/vtk.h"


namespace ale {
namespace mesh {

#ifdef HAVE_VTK 

////////////////////////////////////////////////////////////////////////////////
//! figure out the vtk array type at compile time
////////////////////////////////////////////////////////////////////////////////
template< typename T >
struct vtk_array_t {};

template<>
struct vtk_array_t<float>
{ 
  using type = vtkFloatArray;
};

template<>
struct vtk_array_t<double>
{ 
  using type = vtkDoubleArray;
};

template<>
struct vtk_array_t<int>
{ 
  using type = vtkIntArray;
};

template<>
struct vtk_array_t<long>
{ 
  using type = vtkLongArray;
};

template<>
struct vtk_array_t<long long>
{ 
  using type = vtkLongLongArray;
};

#endif



////////////////////////////////////////////////////////////////////////////////
//! \brief convert a burton mesh to a vtk unstructured grid
////////////////////////////////////////////////////////////////////////////////
template< typename M >
static auto to_vtk( M & m ) 
{

#ifdef HAVE_VTK 

  // alias some types
  using std::vector;

  using   real_t = typename M::real_t;
  using integer_t= typename M::integer_t;
  using vector_t = typename M::vector_t;


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
  auto rspav = access_type_if(m, real_t, is_persistent_at(m,vertices));
  for(auto sf: rspav) {
    auto label = validate_string( sf.label() );      
    auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
    vals->SetNumberOfValues( num_vertices );
    vals->SetName( label.c_str() );
    for(auto v: m.vertices()) vals->SetValue( v.id(), sf[v] );
    pd->AddArray( vals );
  } // for

  // int scalars persistent at vertices
  auto ispav = access_type_if(m, integer_t, is_persistent_at(m,vertices));
  for(auto sf: ispav) {
    auto label = validate_string( sf.label() );
    auto vals = vtkSmartPointer< typename vtk_array_t<integer_t>::type >::New();
    vals->SetNumberOfValues( num_vertices );
    vals->SetName( label.c_str() );
    for(auto v: m.vertices()) vals->SetValue( v.id(), sf[v] );
    pd->AddArray( vals );
  } // for

  // real vectors persistent at vertices

  auto rvpav = access_type_if(m, vector_t, is_persistent_at(m,vertices));
  for(auto vf: rvpav) {
    auto label = validate_string( vf.label() );
    auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
    vals->SetNumberOfComponents( 3 ); // always 3d
    vals->SetNumberOfTuples( num_vertices );
    vals->SetName( label.c_str() );
    for(auto v: m.vertices()) {
      real_t to_vals[3] = {0, 0, 0};
      auto & from_vals = vf[v];
      std::copy( from_vals.begin(), from_vals.end(), to_vals );
      vals->SetTuple( v.id(), to_vals );
    } // for
    pd->AddArray( vals );
  } // for

  //----------------------------------------------------------------------------
  // cell field data header

  // get the point data object
  auto cd = ug->GetCellData();


  // real scalars persistent at cells
  auto rspac = access_type_if(m, real_t, is_persistent_at(m,cells));
  for(auto sf: rspac) {
    auto label = validate_string( sf.label() );      
    auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
    vals->SetNumberOfValues( num_cells );
    vals->SetName( label.c_str() );
    for(auto c: m.cells()) vals->SetValue( c.id(), sf[c] );
    cd->AddArray( vals );
  } // for

  // int scalars persistent at cells
  auto ispac = access_type_if(m, integer_t, is_persistent_at(m,cells));
  for(auto sf: ispac) {
    auto label = validate_string( sf.label() );
    auto vals = vtkSmartPointer< typename vtk_array_t<integer_t>::type >::New();
    vals->SetNumberOfValues( num_cells );
    vals->SetName( label.c_str() );
    for(auto c: m.cells()) vals->SetValue( c.id(), sf[c] );
    cd->AddArray( vals );
  } // for
    
  // real vectors persistent at cells
  auto rvpac = access_type_if(m, vector_t, is_persistent_at(m,cells));
  for(auto vf: rvpac) {
    auto label = validate_string( vf.label() );
    auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
    vals->SetNumberOfComponents( 3 ); // always 3d
    vals->SetNumberOfTuples( num_cells );
    vals->SetName( label.c_str() );
    for(auto c: m.cells()) {
      real_t to_vals[3] = {0, 0, 0};
      auto & from_vals = vf[c];
      std::copy( from_vals.begin(), from_vals.end(), to_vals );
      vals->SetTuple( c.id(), to_vals );
    } // for
    cd->AddArray( vals );
  } // for

  //----------------------------------------------------------------------------
  // return mesh
  //----------------------------------------------------------------------------
  return ug;

#else

    raise_implemented_error( "not built with vtk support." );

    return -1;

#endif

}

////////////////////////////////////////////////////////////////////////////////
//! \brief convert a vtk unstructured grid to a burton mesh
////////////////////////////////////////////////////////////////////////////////
template< typename M >
static auto to_mesh( vtkUnstructuredGrid* ug ) 
{

#ifdef HAVE_VTK 

  // alias some types
  using std::vector;

  using   real_t = typename M::real_t;
  using integer_t= typename M::integer_t;
  using vector_t = typename M::vector_t;
  using  point_t = typename M::point_t;
  using vertex_t = typename M::vertex_t;

  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------
  
  // create a new mesh
  M m;
  
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


  for ( auto v : m.vertices() ) {
    auto & coord = v->coordinates();
  }

  return m;

#else

    raise_implemented_error( "not built with vtk support." );

    return -1;

#endif

}

} // namespace
} // namespace
