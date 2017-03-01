/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file vornoi.h
/// \brief Details for generating vornoit mesh.
/// \remark these are the hidden implementation details
////////////////////////////////////////////////////////////////////////////////
#pragma once


#ifdef HAVE_SHAPO


// library includes
#include <shapo/Mesh2D.hxx>
#include <shapo/Tessellator2D.hxx>
#include <vtkPoints.h>

namespace ale {
namespace mesh {
namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Global type aliases
////////////////////////////////////////////////////////////////////////////////

//! shapo real and index type
using shapo_real_t = double;
using shapo_index_t = shapo::IdType;

////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS - THESE SHOULD NOT BE CALLED OUTSIDE OF THIS UNIT
////////////////////////////////////////////////////////////////////////////////



//==============================================================================
//! \brief a helper function to dump data to a file
//! \param [in] a  the array to dump
//! \param [in] name  the file name to write to
//! \param [in] blksize  the stride of the 2nd array dimension
//==============================================================================
template< class ... Types, 
          template <typename ...> class Container >
static 
void to_file( const Container<Types...> & a, 
              const std::string & name,
              size_t blksize = 1 ) 
{
  std::ofstream myfile( name );

  if (myfile.is_open()) {

    std::size_t cnt = 0;
    std::size_t n = a.size() / blksize;

    myfile << "id = " << 0 << endl;
    myfile << "nb = " << n << endl;
    myfile << "vals = " << endl;

    for ( std::size_t i = 0; i < n; i++ ) {
      for ( std::size_t j = 0; j < blksize; j++ )
        myfile << std::setprecision(17) << std::scientific << a[cnt++] << " ";
      myfile << endl;      
    }
    myfile.close();

  }
  else 
    cerr << "Unable to open file";

}



//==============================================================================
//! \brief create boundaries using all points/edges
//! \param [in] mesh the delfi mesh to extract from
//! \param [in,out] tess the tesselator to inject to
//==============================================================================
template< typename mesh_t >
void extract_boundaries( 
  const mesh_t & mesh, 
  std::vector<shapo_real_t>  & shapo_bnd_point_coords,
  std::vector<shapo_index_t> & shapo_bnd_edge_points ) 
{

  using size_t = typename mesh_t::size_t;
  using point_t = typename mesh_t::point_t;

  // the number of dimenions and max points
  constexpr auto num_dims = mesh_t::num_dimensions;
  auto num_points = mesh.num_vertices();

  // how many boundary points are there
  auto bnd_points = mesh.boundary_vertices();
  auto num_bnd_points = bnd_points.size();

  // create storage for generators
  std::vector<point_t>      bnd_point_coords(num_bnd_points);
  std::map<size_t, size_t>  point_to_bnd_point_id;

  // now set the actual generators 
  auto xp = mesh.vertices()[0]->coordinates();


  // create storage for the min and max boundary points
  auto min_point = xp;
  auto max_point = min_point;

  // now set the actual coordinates
  size_t bid = 0;
  for ( auto v : bnd_points ) {
    // keep track of the index mapping
    point_to_bnd_point_id[v.id()] = bid;
    // get coordinates
    auto xp = v->coordinates();
    // copy coordinates and keep track of min and maxes
    min_point = min( min_point, xp );
    max_point = max( max_point, xp );
    // insert the ponit
    bnd_point_coords[bid++] = std::move(xp);
  }


  // how many boundary edges are there
  auto bnd_edges = mesh.boundary_edges();
  auto num_bnd_edges = bnd_edges.size();

  // create storage for edges
  std::vector< std::pair<size_t,size_t> > bnd_edge_points(num_bnd_edges);

  // now find the actual edges
  bid = 0;
  for ( auto e : bnd_edges ) {
    // Get the two points.  The ordering is consisten because the
    // edge is oriented so the first side is on the inside, and the
    // other is on the outside of the domain.  Ande p1 is p1 for the
    // first side.
    auto points = mesh.vertices(e);
    // the side ordering of points is supposed to be consistent
    bnd_edge_points[bid++] =
      std::make_pair( point_to_bnd_point_id[ points.front().id()  ],
                      point_to_bnd_point_id[ points.back().id() ] );
  }

  // Once the boundary points and edges are extracted, we need to find the 
  // domain corners.
  
  // Just use the bounding box for now.  In the future we will have to
  // walk the edges and find corners ourselves

  constexpr auto num_shapo_bnd_points = 4;

  // Create the list of boundary points
  shapo_bnd_point_coords.reserve( num_dims * num_shapo_bnd_points );

  shapo_bnd_point_coords.emplace_back( min_point[0] );
  shapo_bnd_point_coords.emplace_back( min_point[1] );

  shapo_bnd_point_coords.emplace_back( max_point[0] );
  shapo_bnd_point_coords.emplace_back( min_point[1] );

  shapo_bnd_point_coords.emplace_back( max_point[0] );
  shapo_bnd_point_coords.emplace_back( max_point[1] );

  shapo_bnd_point_coords.emplace_back( min_point[0] );
  shapo_bnd_point_coords.emplace_back( max_point[1] );

  // now create the edges
  shapo_bnd_edge_points = { 
    0, 1, 
    1, 2, 
    2, 3,
    3, 0 
  };

#ifdef DUMP_DIAGNOSTICS
  auto prec = std::cout.precision();
  std::cout.precision(16);
  std::cout.setf ( std::ios::scientific );
  std::cout << " min = (" << min_point[0] << ", " << min_point[1] << ")" << std::endl;
  std::cout << " max = (" << max_point[0] << ", " << max_point[1] << ")" << std::endl;
  std::cout.unsetf ( std::ios::scientific );
  std::cout.precision(prec);
#endif

}

//==============================================================================
//! \brief create the generators
//! \param [in] mesh the delfi mesh to extract from
//! \param [in,out] tess the tesselator to inject to
//==============================================================================
template< typename mesh_t, typename Tessellator >
void create_generators( const mesh_t & mesh, Tessellator & tess ) {

  // the number of dimenions
  constexpr auto num_dims = mesh_t::num_dimensions;

  // figure out how many generators there are:
  // - one generator for each cell
  auto num_gens = mesh.num_cells();

  // create storage for generators
  std::vector<shapo_real_t> generators;
  generators.reserve(num_dims*num_gens);

  // now set the actual generators 
  for( auto c : mesh.cells() ) {
    // get the centroid
    auto xc = c->centroid();
    // push them to the generators
    for ( auto x : xc ) generators.emplace_back( std::move(x) );
  }

  // send the generators to the tessellator
  tess.SetGenerators( generators );

#ifdef DUMP_DIAGNOSTICS
  to_file( generators, "fshapo_setgenerators.dat", num_dims );
#endif

}



//==============================================================================
//! \brief create boundaries using all points/edges
//! \param [in] mesh the delfi mesh to extract from
//! \param [in,out] tess the tesselator to inject to
//==============================================================================
template< typename mesh_t, typename Tessellator >
void create_boundaries_clipped(const mesh_t & mesh, Tessellator & tess) 
{

  // the number of dimenions
  constexpr auto num_dims = mesh_t::num_dimensions;

  // extract the boundaries from the mesh
  std::vector<shapo_real_t>  bnd_point_coords;
  std::vector<shapo_index_t>  bnd_edge_points;
  extract_boundaries( mesh, bnd_point_coords, bnd_edge_points );

  // Send the boundary points to the tessellator
  tess.SetBoundaryPoints( bnd_point_coords );

  // Send the edges to the tessellator
  tess.SetBoundaryEdges( bnd_edge_points );

#ifdef DUMP_DIAGNOSTICS
  to_file( bnd_point_coords, "fshapo_setboundarypoints.dat", num_dims );
  to_file( bnd_edge_points,  "fshapo_setboundaryedges.dat", 2 );
#endif

}



//==============================================================================
//! \brief create boundaries using all points/edges
//! \param [in] mesh the delfi mesh to extract from
//! \param [in,out] tess the tesselator to inject to
//==============================================================================
template< typename mesh_t, typename Tessellator >
void create_boundaries_constrained( const mesh_t & mesh, Tessellator & tess ) 
{

  // the number of dimenions
  constexpr auto num_dims = mesh_t::num_dimensions;

  // the number of generators
  auto num_gens = tess.GetNumberOfGenerators();
  // get the generators
  std::vector<shapo_real_t> generators;
  tess.GetGenerators( generators );

  // extract the boundaries from the mesh
  std::vector<shapo_real_t>  bnd_point_coords;
  std::vector<shapo_index_t>  bnd_edges;
  extract_boundaries( mesh, bnd_point_coords, bnd_edges );

  // add to the list of boundary generators
  for ( auto x : bnd_point_coords )
    generators.emplace_back( std::move(x) );

  // reset the generators
  tess.SetGenerators( generators );

  // build the list of boundary edges.  remeber, we tacked the
  // boundary generators on the end
  for ( auto & i : bnd_edges )
    i += num_gens;
    
  // now set the boundary edges
  tess.SetBoundaryEdges( bnd_edges );

#ifdef DUMP_DIAGNOSTICS
  to_file( generators, "fshapo_setgenerators.dat", num_dims );
  to_file( bnd_edges,  "fshapo_setboundaryedges.dat", 2 );
#endif

}

//==============================================================================
//! \brief transfer points from shapo's mesh to flag's mesh
//! \param [in] tess the tesselator to extract from
//! \param [in] mesh_handle the handle of the parent to create under
//! \return a unique pointer to the mesh
//==============================================================================
template < typename mesh_t, typename Tessellator >
mesh_t transfer_mesh( const Tessellator & tess ) 
{


  // some stl features and usefule types
  using std::swap;
  using std::pair;

  using counter_t= typename mesh_t::counter_t;
  using size_t   = typename mesh_t::size_t;
  using real_t   = typename mesh_t::real_t;
  using point_t  = typename mesh_t::point_t;
  using vertex_t = typename mesh_t::vertex_t;

  // get the shapo mesh
  auto shapo_mesh = tess.GetOutputMesh();
  auto vtk_mesh   = shapo_mesh->GetMesh();

  // the delfi data that needs populated
  mesh_t mesh;

  // the number of dimensions
  constexpr auto num_dims = mesh_t::num_dimensions;

  //----------------------------------------------------------------------------
  // points
  //----------------------------------------------------------------------------


  // how many points
  auto num_points = shapo_mesh->GetNumberOfPoints();

  // create the storage for points
  mesh.init_parameters( num_points );

  // create the individual vertices
  std::vector<vertex_t*> vs( num_points );

  // get the points, always 3d
  auto points  = shapo_mesh->GetPoints();
    
  // store the max and min
  std::array<shapo_real_t, 3> min_point;
  std::array<shapo_real_t, 3> max_point;
  points->GetPoint(0, min_point.data() );
  points->GetPoint(0, max_point.data() );

  // now copy the points
  for( counter_t p=0; p<num_points; ++p ) {
    // get the point
    auto x = points->GetPoint(p);
    // copy to new array
    point_t point;
    for ( int dim = 0; dim < num_dims; dim++ ) {
      point[dim] = x[dim];
      min_point[dim] = std::min( x[dim], min_point[dim] );
      max_point[dim] = std::max( x[dim], max_point[dim] );
    }
    // create a new vertex
    vs[p] = mesh.create_vertex( point );
  }

    
#ifdef DUMP_DIAGNOSTICS
    auto prec = std::cout.precision();
    std::cout.precision(16);
    std::cout.setf ( std::ios::scientific );
    std::cout << " shapo min = (" << min_point[0] << ", " << min_point[1] << ")" << std::endl;
    std::cout << " shapo max = (" << max_point[0] << ", " << max_point[1] << ")" << std::endl;
    std::cout.unsetf ( std::ios::scientific );
    std::cout.precision(prec);
#endif
    


  //----------------------------------------------------------------------------
  // Cells
  //----------------------------------------------------------------------------

  // how many regular cells are there
  auto num_zones  = shapo_mesh->GetNumberOfCells();

  // create the individual vertices
  std::vector<vertex_t*> zone_points;
  zone_points.reserve( num_dims*num_dims );

  // now create the cells
  for( counter_t z=0; z<num_zones; z++ ) {
    // clear the points array
    zone_points.clear();
    // get the number of cells
    auto num_zone_points = shapo_mesh->GetNumberOfCellPoints(z);
    auto zone_point_ids = shapo_mesh->GetCellPoints(z);
    // reserve storage
    zone_points.reserve( num_zone_points );
    // now create the zone point array
    for ( auto p : zone_point_ids )
      zone_points.emplace_back( vs[p] );
    // create the cell
    mesh.create_cell( zone_points );
  }

  // now finalize the mesh setup
  mesh.init();

  return mesh;

} // transfer_mesh

} // namespace detail
} // namespace mesh
} // namespace ale

#endif

