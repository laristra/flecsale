#ifndef FLESCALE_VTKUNSTRUCTURED_GRID
#define FLESCALE_VTKUNSTRUCTURED_GRID


// user includes
#include <flecsale/mesh/mesh_utils.h>
#include <flecsale/utils/time_utils.h>


// system includes
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>


// vtk includes
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkLongLongArray.h>

using namespace apps::hydro;

using mesh_t = mesh_3d_t;
using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;


// Base template from flescale/ale/mesh/vtk_utils.h
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


// a lambda function for validating strings
auto validate_string = []( auto && str ) 
{
	return std::forward<decltype(str)>(str);
};


inline void initVTKUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> grid)
{
	grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
}



template< typename T >
inline void writePointsToGrid(T & mesh, vtkSmartPointer<vtkUnstructuredGrid> grid)
{

	auto points = vtkPoints::New();
	points->SetNumberOfPoints( mesh.num_vertices() );

	for ( auto v : mesh.vertices() ) 
	{
		auto id = v.id();
		auto & coord = v->coordinates();
		double xyz[3] = {0, 0, 0};
		std::copy(coord.begin(), coord.end(), xyz);
		points->SetPoint(id, xyz);
	}

	// transfer the points
	grid->SetPoints(points);
	points->Delete();
}

template< typename T >
inline void writeCellsToGrid(T & mesh, vtkSmartPointer<vtkUnstructuredGrid> grid)
{
	for ( auto c : mesh.cells() ) 
	{
	    // get the vertices in this cell
	    auto cell_verts = mesh.vertices(c);
	    auto num_cell_verts = cell_verts.size();

	    // copy them to the vtk type
	    std::vector< vtkIdType > vert_ids(num_cell_verts);
	    std::transform( cell_verts.begin(), cell_verts.end(), vert_ids.begin(), [](auto && v) { return v.id(); } );

	    // get the faces
	    auto cell_faces = mesh.faces(c);
	    auto num_cell_faces = cell_faces.size();

	    // get the total number of vertices
	    auto tot_verts = std::accumulate( 
	      	cell_faces.begin(), cell_faces.end(), static_cast<size_t>(0),
	      	[&mesh](auto sum, auto f) { return sum + mesh.vertices(f).size(); }
	    );

	    // the list of faces that vtk requires contains the number of points in each
	    // face AND the point ids themselves.
	    std::vector< vtkIdType > face_data;
	    face_data.reserve( tot_verts + num_cell_faces );
	    for ( auto f : cell_faces ) 
	    {
	    	auto face_cells = mesh.cells(f);
	      	auto face_verts = mesh.vertices(f);
	      	auto num_face_verts = face_verts.size();

	      	// copy the face vert ids to the vtk type
	      	std::vector< vtkIdType > face_vert_ids( num_face_verts );
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
	    grid->InsertNextCell(
	      VTK_POLYHEDRON, num_cell_verts, vert_ids.data(),
	      num_cell_faces, face_data.data()
	    );
	}
}



inline void outputVTKUnstructredGrid(std::string filename, vtkSmartPointer<vtkUnstructuredGrid> grid)
{
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInputDataObject(grid);
	writer->SetFileTypeToBinary();
	writer->SetFileName(filename.c_str());
	writer->Write();
}


//inline void populate(auto mesh, vtkSmartPointer<vtkUnstructuredGrid> grid)
template< typename T >
inline void populate(T & mesh, vtkSmartPointer<vtkUnstructuredGrid> grid, unsigned int num_steps)
{
	std::cout << "\n\nWriting to catalyst for ts: " << num_steps << std::endl;

	writePointsToGrid( mesh, grid);
	writeCellsToGrid( mesh, grid);


	
	//
	// 3. Populate fields
	//using   real_t = typename common::real_t;
	//using integer_t= typename common::integer_t;
	//using vector_t = typename mesh_t::vector_t;

	using integer_t = typename mesh_t::integer_t;
  	using real_t = typename mesh_t::real_t;
  	using vector_t = typename mesh_t::vector_t;



	// some general mesh stats
	auto num_vertices = mesh.num_vertices();
	auto num_cells = mesh.num_cells();

	std::cout << "num_vertices: " << num_vertices << std::endl;
	std::cout << "num_cells: " << num_cells << std::endl;

	//----------------------------------------------------------------------------
	// nodal field data header

	// get the point data object
	auto pd = grid->GetPointData();

	
	// real scalars persistent at vertices
	auto rspav = flecsi_get_accessors_all(mesh, real_t, dense, 0, flecsi_has_attribute_at(persistent,vertices));

	for( auto sf: rspav) 
	{
		auto label = validate_string( sf.label() );      
		auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
		vals->SetNumberOfValues( num_vertices );
		vals->SetName( label.c_str() );
		std::cout << "node real scalars: " << label.c_str() << std::endl;
		for(auto v: mesh.vertices()) vals->SetValue( v.id(), sf[v] );
			pd->AddArray( vals );
	}


	// int scalars persistent at vertices
	auto ispav = flecsi_get_accessors_all( mesh, integer_t, dense, 0, flecsi_has_attribute_at(persistent,vertices));

	for(auto sf: ispav) 
	{
		auto label = validate_string( sf.label() );
		auto vals = vtkSmartPointer< typename vtk_array_t<integer_t>::type >::New();
		vals->SetNumberOfValues( num_vertices );
		vals->SetName( label.c_str() );
		std::cout << "node int scalars: " << label.c_str() << std::endl;
		for(auto v: mesh.vertices()) vals->SetValue( v.id(), sf[v] );
			pd->AddArray( vals );
	}

	
	int _count = 0;
	// real vectors persistent at vertices
	auto rvpav = flecsi_get_accessors_all( mesh, vector_t, dense, 0, flecsi_has_attribute_at(persistent,vertices));

	
	for(auto vf: rvpav) 
	{
		std::cout << "_count " << _count << std::endl;
		auto label = validate_string( vf.label() );
		auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
		vals->SetNumberOfComponents( 3 ); // always 3d
		vals->SetNumberOfTuples( num_vertices );
		vals->SetName( label.c_str() );
		std::cout << "node real vectors: " << label.c_str() << std::endl;

		for(auto v: mesh.vertices()) 
		{
			real_t to_vals[3] = {0, 0, 0};
			auto & from_vals = vf[v];
			std::copy( from_vals.begin(), from_vals.end(), to_vals );
			vals->SetTuple( v.id(), to_vals );
		}

		pd->AddArray( vals );
		_count++;
	}

	
	//----------------------------------------------------------------------------
	// cell field data header

	// get the point data object
	auto cd = grid->GetCellData();

	// real scalars persistent at cells
	auto rspac = flecsi_get_accessors_all(mesh, real_t, dense, 0, flecsi_has_attribute_at(persistent,cells));
	for(auto sf: rspac) 
	{
		auto label = validate_string( sf.label() );      
		auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
		vals->SetNumberOfValues( num_cells );
		vals->SetName( label.c_str() );
		std::cout << "cell real scalars: " << label.c_str() << std::endl;
		for(auto c: mesh.cells()) vals->SetValue( c.id(), sf[c] );
			cd->AddArray( vals );
	}


	// int scalars persistent at cells
	auto ispac = flecsi_get_accessors_all(mesh, integer_t, dense, 0, flecsi_has_attribute_at(persistent,cells));
	for(auto sf: ispac) 
	{
		auto label = validate_string( sf.label() );
		auto vals = vtkSmartPointer< typename vtk_array_t<integer_t>::type >::New();
		vals->SetNumberOfValues( num_cells );
		vals->SetName( label.c_str() );
		std::cout << "cell int scalars: " << label.c_str() << std::endl;

		for(auto c: mesh.cells()) vals->SetValue( c.id(), sf[c] );
			cd->AddArray( vals );
	}

	// real vectors persistent at cells
	_count = 0;
	auto rvpac = flecsi_get_accessors_all(mesh, vector_t, dense, 0, flecsi_has_attribute_at(persistent,cells));
	for(auto vf: rvpac) 
	{
		std::cout << "_count " << _count << std::endl;
		auto label = validate_string( vf.label() );
		auto vals = vtkSmartPointer< typename vtk_array_t<real_t>::type >::New();
		vals->SetNumberOfComponents( 3 ); // always 3d
		vals->SetNumberOfTuples( num_cells );
		vals->SetName( label.c_str() );
		std::cout << "cell real vectors: " << label.c_str() << std::endl;
	
		for(auto c: mesh.cells()) 
		{
			real_t to_vals[3] = {0, 0, 0};
			auto & from_vals = vf[c];
			std::copy( from_vals.begin(), from_vals.end(), to_vals );
			vals->SetTuple( c.id(), to_vals );
		} // for

		_count++;
		cd->AddArray( vals );
	}


	//std::string _name = "sep_test__"+ std::to_string(num_steps) +".vtk";
	//outputVTKUnstructredGrid(_name, grid);
}

//} // end of interface
#endif


