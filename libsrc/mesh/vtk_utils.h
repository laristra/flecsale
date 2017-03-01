/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once


// user includes
#include "io/vtk.h"

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


namespace ale {
namespace mesh {

#ifdef HAVE_VTK 

////////////////////////////////////////////////////////////////////////////////
//! \brief A struct to figure out the vtk array type at compile time.
////////////////////////////////////////////////////////////////////////////////
//! @{

//! \remark Base template.
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
//! @}

#else

////////////////////////////////////////////////////////////////////////////////
//! \brief Guess the native vtk integer type if there is no VTK library.
////////////////////////////////////////////////////////////////////////////////
using vtkIdType = int;

#endif


////////////////////////////////////////////////////////////////////////////////
//! \brief The vtk real type
////////////////////////////////////////////////////////////////////////////////
using vtkRealType = double;


#if HAVE_VTK

namespace detail {

////////////////////////////////////////////////////////////////////////////////
//! \brief Write the field data from a mesh to a vtkUnstructuredMesh.
//! \param [in] m A mesh to convert to vtk.
//! \param [in,out] ug A vtk unstructured mesh object.
//! \return a VTK mesh object constructed using \a m.
////////////////////////////////////////////////////////////////////////////////
template< typename M >
static auto write_fields_to_vtk( M & m, vtkUnstructuredGrid* ug ) 
{

  using   real_t = typename M::real_t;
  using integer_t= typename M::integer_t;
  using vector_t = typename M::vector_t;

  // a lambda function for validating strings
  auto validate_string = []( auto && str ) {
    return std::forward<decltype(str)>(str);
  };

  // some general mesh stats
  auto num_vertices = m.num_vertices();
  auto num_cells = m.num_cells();
  
  //----------------------------------------------------------------------------
  // nodal field data header

  // get the point data object
  auto pd = ug->GetPointData();

  // real scalars persistent at vertices
  auto rspav = flecsi_get_accessors_all(
    m, real_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
  );
  for(auto sf: rspav) {
    auto label = validate_string( sf.label() );      
    auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
    vals->SetNumberOfValues( num_vertices );
    vals->SetName( label.c_str() );
    for(auto v: m.vertices()) vals->SetValue( v.id(), sf[v] );
    pd->AddArray( vals );
  } // for

  // int scalars persistent at vertices
  auto ispav = flecsi_get_accessors_all(
    m, integer_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
  );
  for(auto sf: ispav) {
    auto label = validate_string( sf.label() );
    auto vals = vtkSmartPointer< typename vtk_array_t<integer_t>::type >::New();
    vals->SetNumberOfValues( num_vertices );
    vals->SetName( label.c_str() );
    for(auto v: m.vertices()) vals->SetValue( v.id(), sf[v] );
    pd->AddArray( vals );
  } // for

  // real vectors persistent at vertices
  auto rvpav = flecsi_get_accessors_all(
    m, vector_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
  );
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
  auto rspac = flecsi_get_accessors_all(
    m, real_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
  );
  for(auto sf: rspac) {
    auto label = validate_string( sf.label() );      
    auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
    vals->SetNumberOfValues( num_cells );
    vals->SetName( label.c_str() );
    for(auto c: m.cells()) vals->SetValue( c.id(), sf[c] );
    cd->AddArray( vals );
  } // for

  // int scalars persistent at cells
  auto ispac = flecsi_get_accessors_all(
    m, integer_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
  );
  for(auto sf: ispac) {
    auto label = validate_string( sf.label() );
    auto vals = vtkSmartPointer< typename vtk_array_t<integer_t>::type >::New();
    vals->SetNumberOfValues( num_cells );
    vals->SetName( label.c_str() );
    for(auto c: m.cells()) vals->SetValue( c.id(), sf[c] );
    cd->AddArray( vals );
  } // for
    
  // real vectors persistent at cells
  auto rvpac = flecsi_get_accessors_all(
    m, vector_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
  );
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
  return ug;

}


////////////////////////////////////////////////////////////////////////////////
//! \brief Write the vertex coordinates from a mesh to a vtkUnstructuredMesh.
//! \param [in] m A mesh to convert to vtk.
//! \param [in,out] ug A vtk unstructured mesh object.
//! \return a VTK mesh object constructed using \a m.
////////////////////////////////////////////////////////////////////////////////
template< typename M >
static auto write_points_to_vtk( M & m, vtkUnstructuredGrid* ug ) 
{

  using   real_t = typename M::real_t;

  // some general mesh stats
  auto num_vertices = m.num_vertices();

  // figure out the data type
  using data_type = typename vtk_array_t<real_t>::type;
  auto tmp_data = vtkSmartPointer<data_type>::New();

  //----------------------------------------------------------------------------
  // Write the points

  // create points
  auto points = vtkPoints::New();
  points->SetDataType( tmp_data->GetDataType() );
  points->SetNumberOfPoints( num_vertices );
  for ( auto v : m.vertices() ) {
    auto id = v.id();
    auto & coord = v->coordinates();
    vtkRealType x[3] = {0, 0, 0};
    std::copy( coord.begin(), coord.end(), x );
    points->SetPoint( id, x );
  }

  // transfer the points
  ug->SetPoints(points);
  points->Delete();

  //----------------------------------------------------------------------------
  // return mesh
  return ug;

}


} // namespace detail

#endif // HAVE_VTK

////////////////////////////////////////////////////////////////////////////////
//! \brief convert a burton mesh to a vtk unstructured grid
//! \param [in] m A mesh to convert to vtk.
//! \return a VTK mesh object constructed using \a m.
//! \remark This is the 2D version
////////////////////////////////////////////////////////////////////////////////
template< 
  typename M,
  bool enabled = ( M::num_dimensions == 2 ),
  std::enable_if_t< enabled, M >* = nullptr
>
static auto to_vtk( M & m ) 
{

#ifdef HAVE_VTK 

  // alias some types
  using std::vector;

  using   real_t = typename M::real_t;
  using integer_t= typename M::integer_t;
  using vector_t = typename M::vector_t;


  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------

  // some general mesh stats
  auto num_cells = m.num_cells();
    
  // creat unstructured grid
  auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // extract the points
  detail::write_points_to_vtk( m, ug );

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
  detail::write_fields_to_vtk( m, ug );

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
//! \brief convert a burton mesh to a vtk unstructured grid
//! \param [in] m A mesh to convert to vtk.
//! \return a VTK mesh object constructed using \a m.
//! \remark This is the 3D version
////////////////////////////////////////////////////////////////////////////////
template< 
  typename M,
  bool enabled = ( M::num_dimensions == 3 ),
  typename = std::enable_if_t< enabled, M >
>
static auto to_vtk( M & m ) 
{

#ifdef HAVE_VTK 


  // alias some types
  using std::vector;

  using   size_t = typename M::size_t;
  using   real_t = typename M::real_t;
  using integer_t= typename M::integer_t;
  using vector_t = typename M::vector_t;


  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------

  // some general mesh stats
  constexpr auto num_dims = M::num_dimensions;
  auto num_vertices = m.num_vertices();
  auto num_cells = m.num_cells();
    
  // creat unstructured grid
  auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // extract the points
  detail::write_points_to_vtk( m, ug );

  // create the cells
  for ( auto c : m.cells() ) {
    // get the vertices in this cell
    auto cell_verts = m.vertices(c);
    auto num_cell_verts = cell_verts.size();
    // copy them to the vtk type
    vector< vtkIdType > vert_ids(num_cell_verts);
    std::transform( cell_verts.begin(), cell_verts.end(), vert_ids.begin(),
                    [](auto && v) { return v.id(); } );
    // get the faces
    auto cell_faces = m.faces(c);
    auto num_cell_faces = cell_faces.size();
    // get the total number of vertices
    auto tot_verts = std::accumulate( 
      cell_faces.begin(), cell_faces.end(), static_cast<size_t>(0),
      [&m](auto sum, auto f) { return sum + m.vertices(f).size(); }
    );
    // the list of faces that vtk requires contains the number of points in each
    // face AND the point ids themselves.
    vector< vtkIdType > face_data;
    face_data.reserve( tot_verts + num_cell_faces );
    for ( auto f : cell_faces ) {
      auto face_cells = m.cells(f);
      auto face_verts = m.vertices(f);
      auto num_face_verts = face_verts.size();
      // copy the face vert ids to the vtk type
      vector< vtkIdType > face_vert_ids( num_face_verts );
      std::transform( 
        face_verts.begin(), face_verts.end(), face_vert_ids.begin(),
        [](auto && v) { return v.id(); } 
      );
      // check the direction of the vertices
      if ( face_cells[0] != c ) 
        std::reverse( face_vert_ids.begin(), face_vert_ids.end() );
      // now copy them to the global array
      face_data.emplace_back( num_face_verts );
      for ( auto v : face_vert_ids )
        face_data.emplace_back( v );
    }
    // set the cell vertices
    ug->InsertNextCell(
      VTK_POLYHEDRON, num_cell_verts, vert_ids.data(),
      num_cell_faces, face_data.data()
    );
  }
    

  //----------------------------------------------------------------------------
  // write field data
  //----------------------------------------------------------------------------
    
  detail::write_fields_to_vtk( m, ug );


  //----------------------------------------------------------------------------
  // return mesh
  //----------------------------------------------------------------------------
  return ug;


#else

    raise_implemented_error( "not built with vtk support." );

    return -1;

#endif

}

//##############################################################################
//##############################################################################
#ifdef HAVE_VTK 
//##############################################################################
//##############################################################################

////////////////////////////////////////////////////////////////////////////////
//! \brief convert a vtk unstructured grid to a burton mesh
//! \param [in] ug  A mesh in VTK's unstructured format.
//! \return A new mesh of type \a M.
//! \remark This is the 2D version.
////////////////////////////////////////////////////////////////////////////////
template< 
  typename M,
  bool enabled = ( M::num_dimensions == 2 )
>
static 
std::enable_if_t< enabled, M >
to_mesh( vtkUnstructuredGrid* ug ) 
{

  // alias some types
  using std::vector;

  using   real_t = typename M::real_t;
  using integer_t= typename M::integer_t;
  using vector_t = typename M::vector_t;
  using  point_t = typename M::point_t;
  using vertex_t = typename M::vertex_t;
  using counter_t= typename M::counter_t;

  //----------------------------------------------------------------------------
  // setup
  //----------------------------------------------------------------------------
  
  // create a new mesh
  M m;
  
  // some general mesh stats
  constexpr auto num_dims = M::num_dimensions;


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

  for (counter_t i = 0; i < num_vertices; ++i) {
    vtkRealType x[3] = {0, 0, 0};
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

  return m;


}

////////////////////////////////////////////////////////////////////////////////
//! \brief convert a vtk unstructured grid to a burton mesh
//! \param [in] ug  A mesh in VTK's unstructured format.
//! \return A new mesh of type \a M.
//! \remark This is the 3D version.
////////////////////////////////////////////////////////////////////////////////
template< 
  typename M,
  bool enabled = ( M::num_dimensions == 3 ),
  typename = std::enable_if_t< enabled, M >
>
static auto to_mesh( vtkUnstructuredGrid* ug ) 
{

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
  constexpr auto num_dims = M::num_dimensions;

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
    vtkRealType x[3] = {0, 0, 0};
    points->GetPoint( i, x );
    point_t p = { static_cast<real_t>(x[0]), 
                  static_cast<real_t>(x[1]),
                  static_cast<real_t>(x[2]) };
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
  for ( std::size_t i=0; i<num_cells; ++i ) {
    // get the vtk cell data
    vtkIdType num_cell_faces, *cell_vert_ids;
    ug->GetFaceStream( i, num_cell_faces, cell_vert_ids );
    // iterate over each face
    raise_implemented_error("havent finished reading vtk yet");
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

  return m;



}

//##############################################################################
//##############################################################################
#endif // HAVE_VTK
//##############################################################################
//##############################################################################

} // namespace
} // namespace
