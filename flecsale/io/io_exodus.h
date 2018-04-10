/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file 
/// \brief Defines the functionality for the exodus writer and reader.
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <flecsale-config.h>

#include <ristra/assertions/errors.h>

// user includes
#ifdef FLECSALE_ENABLE_EXODUS
#  include <flecsi-sp/io/exodus_definition.h>
#endif

// Paraview has a problem with regions in nfaced data.  Uncomment the next
// line, or compile with -dPARAVIEW_EXODUS_3D_REGION_BUGFIX to outout exodus
// files with only one region.
// #define PARAVIEW_EXODUS_3D_REGION_BUGFIX


namespace flecsale {
namespace io {


////////////////////////////////////////////////////////////////////////////////
/// \brief provides base functionality for exodus writer
/// \tparam N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TYPE >
class io_exodus__ {

public:

  //============================================================================
  //! Typedefs
  //============================================================================

  //! the mesh type
  using mesh_t = MESH_TYPE;

  // other useful types
  using    size_t = typename mesh_t::size_t;
  using counter_t = typename mesh_t::counter_t;
  using integer_t = typename mesh_t::integer_t;
  using    real_t = typename mesh_t::real_t;
  using   point_t = typename mesh_t::point_t;
  using  vector_t = typename mesh_t::vector_t;
  using  vertex_t = typename mesh_t::vertex_t;
  using    face_t = typename mesh_t::face_t;
  using    cell_t = typename mesh_t::cell_t;
  using  ex_real_t  = real_t;
  using  ex_index_t = int;


#ifdef FLECSALE_ENABLE_EXODUS

  //! use flecsi's base functionality
  using base_t =
    flecsi_sp::io::exodus_base__< MESH_TYPE::num_dimensions, real_t >;

  //============================================================================
  //! \brief write field data to the file
  //! \param [in] m  The mesh to extract field data from.
  //! \return the status of the file
  //============================================================================
  template< typename T >
  static
  void write_fields( int exoid, mesh_t & m, size_t time_step, const T & f ) 
  { 

    int status;

    // mesh statistics
    constexpr auto num_dims  = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();
    auto num_elem = m.num_cells();
    auto num_elem_blk = m.num_regions();

    //--------------------------------------------------------------------------
    // initial setup

    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    //--------------------------------------------------------------------------
    // hard coded data
    //--------------------------------------------------------------------------

    int num_var = 1;
    int var_id = 1;
    int elem_blk_id = 1;

    // put the number of element fields
    status = ex_put_var_param(exoid, "e", num_var);
    if (status)
      throw_runtime_error(
        "Problem writing variable number, " <<
        " ex_put_var_param() returned " << status 
      );

    // fill element variable names array
    //doesnt work
    //auto label = validate_string( f.label() );
    auto label = std::string( "output_variable" );
    status = ex_put_var_name(exoid, "e", var_id, "density");
    if (status)
      throw_runtime_error(
        "Problem writing variable name, " <<
        " ex_put_var_name() returned " << status 
      );

    size_t cid = 0;
    std::vector<ex_real_t> tmp(m.num_cells()); 
    for(auto c: m.cells()) tmp[cid++] = f(c);
    status = ex_put_elem_var(
      exoid, time_step, var_id, elem_blk_id, tmp.size(), tmp.data()
    );
    if (status)
      throw_runtime_error(
        "Problem writing variable data, " <<
        " ex_put_elem_var() returned " << status 
      );

#if 0

    //--------------------------------------------------------------------------
    // nodal data
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // nodal field data header        

    // number of nodal fields
    int num_nf = 0;

    // real scalars persistent at vertices
    auto rspav = flecsi_get_accessors_all(
      m, real_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
    );
    num_nf += rspav.size();
    // int scalars persistent at vertices
    auto ispav = flecsi_get_accessors_all(
      m, integer_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
    );
    num_nf += ispav.size();
    // real vectors persistent at vertices
    auto rvpav = flecsi_get_accessors_all(
      m, vector_t, dense, 0, flecsi_has_attribute_at(persistent,vertices)
    );
    num_nf += num_dims*rvpav.size();

    // variable extension for vectors
    std::string var_ext[3];
    var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

    // put the number of nodal fields
    if (num_nf > 0) {
      status = ex_put_var_param(exoid_, "n", num_nf);
      assert(status == 0);
    }

    // fill node variable names array
    int inum = 1;

    for(auto sf: rspav) {
      auto label = validate_string( sf.label() );      
      status = ex_put_var_name(exoid_, "n", inum++, label.c_str());
      assert(status == 0);
    } // for
    for(auto sf: ispav) {
      auto label = validate_string( sf.label() );
      status = ex_put_var_name(exoid_, "n", inum++, label.c_str());
    } // for
    for(auto vf: rvpav) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        status = ex_put_var_name(exoid_, "n", inum++, dim_label.c_str());
      } // for
    } // for


    //--------------------------------------------------------------------------
    // nodal field data

    std::vector<ex_real_t> tmp(num_nodes);

    // write nodal fields
    inum = 1;

    // node field buffer
    for(auto sf: rspav) {
      for(auto v: m.vertices()) tmp[v.id()] = sf[v];
      status = ex_put_nodal_var(exoid_, time_step, inum++, num_nodes, tmp.data());
      assert(status == 0);
    } // for
    for(auto sf: ispav) {
      // cast int fields to real_t
      for(auto v: m.vertices()) tmp[v.id()] = (real_t)sf[v];
      status = ex_put_nodal_var(exoid_, time_step, inum++, num_nodes, tmp.data());
      assert(status == 0);
    } // for
    for(auto vf: rvpav) {
      for(int d=0; d < num_dims; ++d) {
        for(auto v: m.vertices()) tmp[v.id()] = vf[v][d];
        status = ex_put_nodal_var(exoid_, time_step, inum++, num_nodes, tmp.data());
        assert(status == 0);
      } // for
    } // for


    //--------------------------------------------------------------------------
    // element fields
    //--------------------------------------------------------------------------

    // get the regions
    auto region_cells = m.regions();

    //--------------------------------------------------------------------------
    // element field data headers

    // number of element fields
    int num_ef = 0;

    // real scalars persistent at cells
    auto rspac = flecsi_get_accessors_all(
      m, real_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
    );
    num_ef += rspac.size();
    // int scalars persistent at cells
    auto ispac = flecsi_get_accessors_all(
      m, integer_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
    );
    num_ef += ispac.size();
    // real vectors persistent at cells
    auto rvpac = flecsi_get_accessors_all(
      m, vector_t, dense, 0, flecsi_has_attribute_at(persistent,cells)
    );
    num_ef += num_dims*rvpac.size();

    // put the number of element fields
    if(num_ef > 0) {
      status = ex_put_var_param(exoid_, "e", num_ef);
      assert(status == 0);
    } // if

    // fill element variable names array
    inum = 1;

    for(auto sf: rspac) {
      auto label = validate_string( sf.label() );
      status = ex_put_var_name(exoid_, "e", inum++, label.c_str());
      assert(status == 0);
    } // for
    for(auto sf: ispac) {
      auto label = validate_string( sf.label() );
      status = ex_put_var_name(exoid_, "e", inum++, label.c_str());
      assert(status == 0);
    } // for
    for(auto vf: rvpac) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        status = ex_put_var_name(exoid_, "e", inum++, dim_label.c_str());
      } // for
    } // for


    //--------------------------------------------------------------------------
    // element field data

    // loop over element blocks, writing the filtered element data
    for ( int iblk=0; iblk<num_elem_blk; iblk++ ) {

      // stats for this block
      auto elem_blk_id = iblk+1;

      // get the elements in this block
      const auto & elem_this_blk = region_cells[iblk];
      auto num_elem_this_blk = elem_this_blk.size();
  
      // resize temp data
      tmp.clear();
      tmp.resize( num_elem_this_blk );

      // write element fields
      inum = 1;

      // element field buffer
      for(auto sf: rspac) {
        size_t cid = 0;
        for(auto c: elem_this_blk) tmp[cid++] = sf[c];
        status = ex_put_elem_var(exoid_, time_step, inum++, elem_blk_id, num_elem_this_blk, tmp.data());
        assert(status == 0);
      } // for
      for(auto sf: ispac) {
        // cast int fields to real_t
        size_t cid = 0;
        for(auto c: elem_this_blk) tmp[cid++] = (real_t)sf[c];
        status = ex_put_elem_var(exoid_, time_step, inum++, elem_blk_id, num_elem_this_blk, tmp.data());
        assert(status == 0);
      } // for
      for(auto vf: rvpac) {
        for(int d=0; d < num_dims; ++d) {
          size_t cid = 0;
          for(auto c: elem_this_blk) tmp[cid++] = vf[c][d];
          status = ex_put_elem_var(exoid_, time_step, inum++, elem_blk_id, num_elem_this_blk, tmp.data());
          assert(status == 0);
        } // for
      } // for

    } // block

    //--------------------------------------------------------------------------
    // final setup

    // set the time
    ex_real_t soln_time = m.time();
    status = ex_put_time(exoid_, time_step, &soln_time );
    assert(status == 0);

#endif 

  }

#endif


  //============================================================================
  //!  \brief Implementation of exodus mesh write for burton specialization.
  //!
  //!  \param[in] name Write burton mesh \e m to \e name.
  //!  \param[in] m Burton mesh to write to \e name.
  //!
  //!  \return Exodus error code. 0 on success.
  //============================================================================

  template< typename T >
  static void write(
    const std::string &name,
    mesh_t &m,
    size_t iteration = 0,
    ex_real_t time = 0.0,
    T * const d = nullptr
  ) {

#ifdef FLECSALE_ENABLE_EXODUS

    std::cout << "Writing mesh to: " << name << std::endl;

    // alias some types
    using std::vector;
    using std::array;


    //--------------------------------------------------------------------------
    // initial setup
    //--------------------------------------------------------------------------
    
    // get the iteration number
    // - first step starts at 1
    // - ignore input iteration because we are always only outputting one solution
    iteration = 1;

    auto exoid = base_t::open( name, std::ios_base::out );
    assert(exoid >= 0);

    // get the general statistics
    constexpr auto num_dims = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();
    auto num_faces = num_dims==3 ? m.num_faces() : 0;
    auto num_elems = m.num_cells();
    auto num_elem_blk = 1;
    auto num_node_sets = 0;
    auto num_side_sets = 0;

    auto exo_params = base_t::make_params();
    exo_params.num_nodes = num_nodes;
    exo_params.num_face = num_faces;
    exo_params.num_face_blk = num_dims==3 ? 1 : 0;;
    exo_params.num_elem = num_elems;
    exo_params.num_elem_blk = 1;
    exo_params.num_node_sets = 0;

    base_t::write_params(exoid, exo_params);

    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    //--------------------------------------------------------------------------
    // Point Coordinates
    //--------------------------------------------------------------------------
    
    std::vector<real_t> vertex_coord( num_nodes * num_dims );
    
    for (auto v : m.vertices()) {
      auto & coords = v->coordinates();
      for ( int i=0; i<num_dims; i++ )
        vertex_coord[ i*num_nodes + v.id() ] = coords[i];
    } // for

    base_t::write_point_coords( exoid, vertex_coord );

  
    //--------------------------------------------------------------------------
    // Face connectivity
    //--------------------------------------------------------------------------

    if ( num_dims == 2 ) {

      // get the master entity lists
      const auto & cs = m.cells();

      // create the element blocks
      base_t::template write_element_block<ex_index_t>( 
        exoid, 1, "cells", num_elems,
        [&]( auto c, auto & face_list ) {
          for ( auto v : m.vertices(cs[c]) ) 
            face_list.emplace_back( v.id() );
        }
      );

    }

    else if ( num_dims == 3 ) {
      
      // get the master entity lists
      const auto & fs = m.faces();
      const auto & cs = m.cells();

      // create the face blocks
      base_t::template write_face_block<ex_index_t>( 
        exoid, 1, "faces", num_faces,
        [&]( auto f, auto & face_conn ) 
        {
          for ( auto v : m.vertices(fs[f]) )
            face_conn.emplace_back( v.id() );
        }
      );


      // create the element blocks
      base_t::template write_element_block<ex_index_t>( 
        exoid, 1, "cells", num_elems,
        [&]( auto c, auto & face_list ) {
          for ( auto f : m.faces(cs[c]) ) 
            face_list.emplace_back( f.id() );
        }
      );

    }

    //--------------------------------------------------------------------------
    // write field data
    //--------------------------------------------------------------------------
    if ( d ) 
      write_fields( exoid, m, iteration, *d );


    //--------------------------------------------------------------------------
    // final setup
    //--------------------------------------------------------------------------

    // set the time
    auto status = ex_put_time(exoid, iteration, &time );
    assert(status == 0);

    //--------------------------------------------------------------------------
    // final step
    //--------------------------------------------------------------------------
    base_t::close( exoid );


#else

    std::cerr << "FLECSI not build with exodus support." << std::endl;

#endif

  } // io_exodus_t::write


}; // struct io_exodus_t


} // namespace io
} // namespace flecsale
