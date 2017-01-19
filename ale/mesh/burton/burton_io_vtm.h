/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Provides functionality for writing vtm files.
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
#  include <vtkCompositeDataIterator.h>
#  include <vtkMultiBlockDataSet.h>
#  include <vtkXMLMultiBlockDataReader.h>
#  include <vtkXMLMultiBlockDataWriter.h>
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
/// \brief This is the mesh reader and writer based on the vtm format.
////////////////////////////////////////////////////////////////////////////////
class burton_io_vtm_t : public flecsi::io::io_base_t<burton_mesh_2d_t> {

public:

  //! Default constructor
  burton_io_vtm_t() {}

  //============================================================================
  //! \brief Implementation of vtm mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return vtm error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_2d_t &m ) override
  {

#ifdef HAVE_VTK

    std::cout << "Writing mesh to: " << name << std::endl;

    // alias some types
    using std::string;
    using std::vector;

    using   mesh_t = burton_mesh_2d_t;
    using   size_t = typename mesh_t::size_t;
    using   real_t = typename mesh_t::real_t;
    using integer_t= typename mesh_t::integer_t;
    using vector_t = typename mesh_t::vector_t;
    using vertex_t = typename mesh_t::vertex_t;


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    //--------------------------------------------------------------------------
    // Write Blocks
    //--------------------------------------------------------------------------

    // get the regions
    auto region_cells = m.regions();
    auto num_regions = m.num_regions();

    // vtkMultiBlockDataSet respresents multi-block datasets. See
    // the class documentation for more information.
    auto mb = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    
    for ( int iblk=0; iblk<num_regions; iblk++ ) {

      // creat unstructured grid
      auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

      //------------------------------------------------------------------------
      // Connectivity

      // the cells for this bloci
      const auto & cells_this_block = region_cells[iblk];
      auto num_cells_this_block = cells_this_block.size();

      // get points for this block
      std::unordered_map< vertex_t*, size_t > points_this_block;
      for ( auto c : cells_this_block )
        for ( auto v : m.vertices(c) ) {
          auto id = points_this_block.size();;
          points_this_block.emplace( std::make_pair( v, id ) );
        }
      auto num_vertices = points_this_block.size(); 

      // create point storage
      auto points = vtkPoints::New();
      points->SetNumberOfPoints( num_vertices );
      
      // now create the vtk points, and create local ids
      for ( auto v : points_this_block ) {
        auto id =  v.second;
        auto coord = v.first->coordinates();
        real_t x[3] = {0, 0, 0};
        std::copy( coord.begin(), coord.end(), x );
        points->SetPoint( id, x );
      }

      // transfer the points
      ug->SetPoints(points);
      points->Delete();

      // create the cells
      for ( auto c : cells_this_block ) {
        // get the vertices in this cell
        auto vs = m.vertices(c);
        auto n = vs.size();
        // copy them to the vtk type
        std::vector< vtkIdType > ids(n);
        std::transform( vs.begin(), vs.end(), ids.begin(),
          [&](auto && v) { return points_this_block.at(v); } );
        // set the cell vertices
        ug->InsertNextCell(VTK_POLYGON, n, ids.data());
      }

      //------------------------------------------------------------------------
      // nodal field data

      // get the point data object
      auto pd = ug->GetPointData();

      // real scalars persistent at vertices
      auto rvals = vtk_array_t<real_t>::type::New();
      rvals->SetNumberOfValues( num_vertices );
      
      auto rspav = get_accessors_all(
        m, real_t, dense, 0, has_attribute_at(persistent,vertices)
      );
      for(auto sf: rspav) {
        auto label = validate_string( sf.label() );      
        rvals->SetName( label.c_str() );
        auto vid = 0;
        for(auto v : points_this_block) rvals->SetValue( vid++, sf[v.first] );
        pd->AddArray( rvals );
      } // for

      rvals->Delete();

      // int scalars persistent at vertices
      auto ivals = vtk_array_t<integer_t>::type::New();
      ivals->SetNumberOfValues( num_vertices );

      auto ispav = get_accessors_all(
        m, integer_t, dense, 0, has_attribute_at(persistent,vertices)
      );
      for(auto sf: ispav) {
        auto label = validate_string( sf.label() );
        ivals->SetName( label.c_str() );
        auto vid = 0;
        for(auto v: points_this_block) ivals->SetValue( vid++, sf[v.first] );
        pd->AddArray( ivals );
      } // for
    
      ivals->Delete();

      // real vectors persistent at vertices
      auto vvals = vtk_array_t<real_t>::type::New();
      vvals->SetNumberOfComponents( 3 ); // always 3d
      vvals->SetNumberOfTuples( num_vertices );

      auto rvpav = get_accessors_all(
        m, vector_t, dense, 0, has_attribute_at(persistent,vertices)
      );
      for(auto vf: rvpav) {
        auto label = validate_string( vf.label() );
        vvals->SetName( label.c_str() );
        auto vid = 0;
        for(auto v: points_this_block) {
          real_t to_vals[3] = {0, 0, 0};
          auto & from_vals = vf[v.first];
          std::copy( from_vals.begin(), from_vals.end(), to_vals );
          vvals->SetTuple( vid++, to_vals );
        } // for
        pd->AddArray( vvals );
      } // for

      vvals->Delete();


      //------------------------------------------------------------------------
      // cell field data header

      // get the point data object
      auto cd = ug->GetCellData();

      // real scalars persistent at cells
      rvals = vtk_array_t<real_t>::type::New();
      rvals->SetNumberOfValues( num_cells_this_block );

      auto rspac = get_accessors_all(
        m, real_t, dense, 0, has_attribute_at(persistent,cells)
      );
      for(auto sf: rspac) {
        auto label = validate_string( sf.label() );      
        rvals->SetName( label.c_str() );
        size_t cid = 0;
        for(auto c: cells_this_block) rvals->SetValue( cid++, sf[c] );
        cd->AddArray( rvals );
      } // for

      rvals->Delete();

      // int scalars persistent at cells
      ivals = vtk_array_t<integer_t>::type::New();
      ivals->SetNumberOfValues( num_cells_this_block );

      auto ispac = get_accessors_all(
        m, integer_t, dense, 0, has_attribute_at(persistent,cells)
      );
      for(auto sf: ispac) {
        auto label = validate_string( sf.label() );
        ivals->SetName( label.c_str() );
        size_t cid = 0;
        for(auto c: cells_this_block) ivals->SetValue( cid++, sf[c] );
        cd->AddArray( ivals );
      } // for
    
      ivals->Delete();

      // real vectors persistent at cells
      vvals = vtk_array_t<real_t>::type::New();
      vvals->SetNumberOfComponents( 3 ); // always 3d
      vvals->SetNumberOfTuples( num_cells_this_block );

      auto rvpac = get_accessors_all(
        m, vector_t, dense, 0, has_attribute_at(persistent,cells)
      );
      for(auto vf: rvpac) {
        auto label = validate_string( vf.label() );
        vvals->SetName( label.c_str() );
        size_t cid = 0;
        for(auto c: cells_this_block) {
          real_t to_vals[3] = {0, 0, 0};
          auto & from_vals = vf[c];
          std::copy( from_vals.begin(), from_vals.end(), to_vals );
          vvals->SetTuple( cid++, to_vals );
        } // for
        cd->AddArray( vvals );
      } // for

      vvals->Delete();

      //------------------------------------------------------------------------
      // Finalize block

      // Add the structured grid to the multi-block dataset
      mb->SetBlock(iblk, ug);
    
    } // block

    // write to file
    auto writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
    writer->SetFileName( name.c_str() );
    writer->SetInputDataObject( mb );

    writer->Write();

    return 0;

#else

    std::cerr << "FLECSI not build with vtk support." << std::endl;
    std::exit(1);

    return -1;

#endif


  }

  //============================================================================
  //! \brief Implementation of vtu mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return vtu error code. 0 on success.
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_2d_t &m) override
  {
#ifdef HAVE_VTK

    std::cout << "Reading mesh from: " << name << std::endl;

    // alias some types
    using std::string;
    using std::vector;

    using   mesh_t = burton_mesh_2d_t;
    using   real_t = typename mesh_t::real_t;
    using integer_t= typename mesh_t::integer_t;
    using vector_t = typename mesh_t::vector_t;
    using  point_t = typename mesh_t::point_t;
    using vertex_t = typename mesh_t::vertex_t;
    using counter_t= typename mesh_t::counter_t;

    constexpr auto test_tolerance = common::test_tolerance;

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    // some general mesh stats
    auto num_dims = m.num_dimensions;

    // a points comparison function
    auto is_same_point = [=](auto first1, auto last1, auto first2) 
      { 
        std::decay_t< decltype(*first1) > dist_sqr(0);
        dist_sqr = std::inner_product( first1, last1, first2, dist_sqr,
          [](auto a, auto b)  { return a + b; },
          [](auto a, auto b)  
          { 
            auto delta = a - b;
            return delta*delta; 
          }
        );
        return ( dist_sqr < test_tolerance );
      };

    //--------------------------------------------------------------------------
    // Read solution
    //--------------------------------------------------------------------------


    // write to file
    auto reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    reader->SetFileName( name.c_str() );
    reader->Update();
    auto mb = reader->GetOutput();


    // get the number of points
    auto num_points = mb->GetNumberOfPoints();

    // storage for vertices
    vector<point_t> ps;
    ps.reserve( num_points );
      
    //--------------------------------------------------------------------------
    // Loop over blocks - Create points
    //--------------------------------------------------------------------------

    // points will overlap at the block boundaries

    auto bit = mb->NewIterator();
    bit->InitTraversal();
    
    auto num_blocks = 0;
    while( !bit->IsDoneWithTraversal() ) {
      
      // get the next block
      auto ug = vtkUnstructuredGrid::SafeDownCast( mb->GetDataSet( bit ) );
      bit->GoToNextItem();
    
      // get points
      auto points_this_block = ug->GetPoints();
      auto num_points_this_block = points_this_block->GetNumberOfPoints();

      // create points
      for (counter_t i = 0; i < num_points_this_block; ++i) {
        vtkRealType x[3] = {0, 0, 0};
        points_this_block->GetPoint( i, x );
        point_t p = { static_cast<real_t>(x[0]), 
                      static_cast<real_t>(x[1]) };
        ps.emplace_back( std::move(p) );
      } // for
            
      // increment block counter
      num_blocks++;

    }  // blocks

    //--------------------------------------------------------------------------
    // Cull matching points
    //--------------------------------------------------------------------------


    // a lexical points comparison function
    auto sort_points = [&num_dims](auto & a, auto & b) 
      { 
        for ( auto d=0; d<num_dims-1; d++ ) 
          if ( a[d] != b[d] ) return ( a[d] < b[d] );
        return ( a[num_dims-1] < b[num_dims-1] );
      };


    // sort the points
    std::sort( ps.begin(), ps.end(), sort_points );

    // a points points comparison function
    auto compare_points = [=](auto & a, auto & b) 
      { return is_same_point(a.begin(), a.end(), b.begin()); };

    // cull the duplicate points
    auto ps_end = std::unique( ps.begin(), ps.end(), compare_points );
    ps.erase( ps_end, ps.end() );

    // tell mesh how many points
    auto num_unique_points = ps.size();
    m.init_parameters( num_unique_points );
    
    // storage for vertices
    vector<vertex_t*> vs;
    vs.reserve( num_unique_points );

    // create the mesh vertices
    for (counter_t i = 0; i < num_unique_points; ++i) {
      auto v = m.create_vertex( ps[i] );
      vs.emplace_back( std::move(v) );
    } // for



    //--------------------------------------------------------------------------
    // Loop over blocks - Create cells
    //--------------------------------------------------------------------------

    // block id counter
    auto block_id = 0;

    // storage for regions
    vector<size_t> region_ids;    
    

    bit->InitTraversal();    
    while( !bit->IsDoneWithTraversal() ) {
      
      // get the next block
      auto ug = vtkUnstructuredGrid::SafeDownCast( mb->GetDataSet( bit ) );
      bit->GoToNextItem();

      // get block points
      auto points_this_block = ug->GetPoints();
      auto num_points_this_block = points_this_block->GetNumberOfPoints();

      // get block cells
      auto cells_this_block = ug->GetCells();
      auto num_cells_this_block = cells_this_block->GetNumberOfCells();

      // storage for element vertices
      vector<vertex_t *> elem_vs;
      vtkIdType npts, *pts;
      
      //------------------------------------------------------------------------
      // loop over cells

      cells_this_block->InitTraversal();
      while( cells_this_block->GetNextCell(npts, pts) ) {
        // reset storage
        elem_vs.clear();
        elem_vs.reserve( npts );
        // loop over each local point
        for ( counter_t v=0;  v<npts; v++ ) {
          //--------------------------------------------------------------------
          // find the closest point
          auto it = std::find_if( vs.begin(), vs.end(),
            [&](auto & a) 
            {
              auto test_coord = a->coordinates();
              vtkRealType coord[3] = {0, 0, 0}; // always 3d
              points_this_block->GetPoint( pts[v], coord );
              return is_same_point(test_coord.begin(), test_coord.end(), coord);
            }
          );
          //--------------------------------------------------------------------
          assert( it != vs.end() && "couldn't find a vertex" ) ;
          elem_vs.emplace_back( *it );
        }
        // create acual cell
        auto c = m.create_cell( elem_vs );                
      }
      // end cells
      //------------------------------------------------------------------------


      // set element regions
      region_ids.resize( region_ids.size() + num_cells_this_block, block_id );

      
      // increment block id counter
      block_id++;

    }  // blocks


    //--------------------------------------------------------------------------
    // Finalize mesh
    //--------------------------------------------------------------------------

    m.init();

    // override the region ids
    for ( auto c : m.cells() )
      c->region() = region_ids[c.id()];
    m.set_num_regions( num_blocks );

    //--------------------------------------------------------------------------
    // Clean up
    //--------------------------------------------------------------------------

    bit->Delete();

    return 0;

#else

    std::cerr << "FLECSI not build with vtk support." << std::endl;
    std::exit(1);

    return -1;

#endif

  };

}; // struct io_vtm_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_vtm_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_vtm_t.
//!
//! \return Pointer to io_base_t base class of io_vtm_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io::io_base_t<burton_mesh_2d_t> * create_io_vtm()
{
  return new burton_io_vtm_t;
} // create_io_vtm


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "vtm" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_vtm_registered =
  flecsi::io::io_factory_t<burton_mesh_2d_t>::instance().registerType(
    "vtm", create_io_vtm );


} // namespace mesh
} // namespace ale
