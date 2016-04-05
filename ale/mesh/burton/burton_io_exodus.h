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
 * \file io_exodus.h
 * \date Initial file creation: Oct 07, 2015
 *
 ******************************************************************************/

#pragma once

//! system includes
#include <cstring>

#ifdef HAVE_EXODUS
#  include <exodusII.h>
#endif

//! user includes
#include "flecsi/io/io_base.h"
#include "../../mesh/burton/burton_mesh.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
/// \class io_exodus_t io_exodus.h
/// \brief io_exodus_t provides a derived type of io_base.h and registrations
///   of the exodus file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_exodus_t : public flecsi::io_base_t<burton_mesh_t> {

  //! Default constructor
  burton_io_exodus_t() {}

  //============================================================================
  //! Implementation of exodus mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_t &m) 
  {

#ifdef HAVE_EXODUS

    std::cout << "Reading mesh from: " << name << std::endl;

    // alias some types
    using std::vector;
    using   mesh_t = burton_mesh_t;
    using   real_t = typename mesh_t::real_t;
    using  point_t = typename mesh_t::point_t;
    using vector_t = typename mesh_t::vector_t;
    using vertex_t = typename mesh_t::vertex_t;
    using ex_real_t = real_t;
    using ex_index_t = int;

    // size of floating point variables used in app.
    int CPU_word_size = sizeof(real_t);
    // size of floating point stored in name.
    int IO_word_size = 0;
    float version;
    
    auto exoid = ex_open(
      name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version);
    assert(exoid >= 0);

    // get the initialization parameters
    ex_init_params exopar;
    auto status = ex_get_init_ext(exoid, &exopar);
    assert(status == 0);

    // verify mesh dimension
    assert(m.num_dimensions() == exopar.num_dim);
    auto num_dims = exopar.num_dim;
    auto num_nodes = exopar.num_nodes;
    auto num_elem = exopar.num_elem;
    auto num_elem_blk = exopar.num_elem_blk;
    auto num_node_sets = exopar.num_node_sets;
    auto num_side_sets = exopar.num_side_sets;

    m.init_parameters(num_nodes);

    // read nodes
    vector<ex_real_t> xcoord(num_nodes);
    vector<ex_real_t> ycoord(num_nodes);
    status = ex_get_coord(exoid, xcoord.data(), ycoord.data(), nullptr);
    assert(status == 0);

    // put nodes into mesh
    vector<vertex_t *> vs;
    vs.reserve( num_nodes );
    for (size_t i = 0; i < num_nodes; ++i) {
      point_t p = { static_cast<real_t>(xcoord[i]), 
                    static_cast<real_t>(ycoord[i]) };
      auto v = m.create_vertex(p);
      vs.push_back(v);
    } // for

    // read blocks
    vector<ex_index_t> blockids( num_elem_blk );
    status = ex_get_elem_blk_ids(exoid, blockids.data() );
    assert(status == 0);

    // storage for regions
    vector<ex_index_t> region_ids;
    region_ids.reserve( num_elem );

    //--------------------------------------------------------------------------
    // read each block
    for ( auto iblk=0; iblk<exopar.num_elem_blk; iblk++ ) {

      auto elem_blk_id = blockids[iblk];

      char block_name[MAX_LINE_LENGTH];
      status = ex_get_name(exoid, EX_ELEM_BLOCK, elem_blk_id, block_name);
      assert(status == 0);

      // get the info about this block
      ex_index_t num_elem_this_blk = 0;
      ex_index_t num_attr = 0;
      ex_index_t num_nodes_per_elem = 0;
      char elem_type[MAX_STR_LENGTH];
      status = ex_get_elem_block(
        exoid, elem_blk_id, elem_type, &num_elem_this_blk, &num_nodes_per_elem, &num_attr);
      assert(status == 0);
      
      //--------------------------------
      // polygon data
      if ( strcmp("nsided",elem_type) == 0 || strcmp("NSIDED",elem_type) == 0 ) {

        // the number of nodes per element is really the number of nodes in the whole block
        auto num_nodes_this_blk = num_nodes_per_elem;


        // get the number of nodes per element
        vector<ex_index_t> elem_node_counts(num_elem_this_blk);
        status = ex_get_entity_count_per_polyhedra( 
          exoid, EX_ELEM_BLOCK, elem_blk_id, elem_node_counts.data() );
        assert(status == 0);

        // read element definitions
        vector<ex_index_t> elem_nodes(num_nodes_this_blk);
        status = ex_get_elem_conn(exoid, elem_blk_id, elem_nodes.data());
        assert(status == 0);
        
        // storage for element verts
        vector<vertex_t *> elem_vs;
        elem_vs.reserve( num_dims * elem_node_counts[0] );
        
        // create cells in mesh
        auto base = 0;
        for (size_t e = 0; e < num_elem_this_blk; ++e) {
          elem_vs.clear();
          // get the number of nodes
          num_nodes_per_elem = elem_node_counts[e];
          // copy local vertices into vector ( exodus uses 1 indexed arrays )
          for ( auto v=0;  v<num_nodes_per_elem; v++ )
            elem_vs.emplace_back( vs[ elem_nodes[base+v] - 1 ] );
          // create acual cell
          auto c = m.create_cell( elem_vs );
          // base offset into elt_conn
          base += num_nodes_per_elem;
        }

      }
      //--------------------------------
      // fixed element size
      else {

        // read element definitions
        vector<ex_index_t> elt_conn(num_elem_this_blk * num_nodes_per_elem);
        status = ex_get_elem_conn(exoid, elem_blk_id, elt_conn.data());
        assert(status == 0);
        
        // storage for element verts
        vector<vertex_t *> elem_vs;
        elem_vs.reserve( num_nodes_per_elem );
        
        // create cells in mesh
        for (size_t e = 0; e < num_elem_this_blk; ++e) {
          elem_vs.clear();
          // base offset into elt_conn
          auto b = e*num_nodes_per_elem;
          // copy local vertices into vector ( exodus uses 1 indexed arrays )
          for ( auto v=0;  v<num_nodes_per_elem; v++ )
            elem_vs.emplace_back( vs[ elt_conn[b+v]-1 ] );
          // create acual cell
          auto c = m.create_cell( elem_vs );          
        }

      } // element type
      //--------------------------------


      // set element regions
      for ( auto e = 0; e < num_elem_this_blk; e++ )
        region_ids.emplace_back( iblk );
      
    }
    // end blocks
    //--------------------------------------------------------------------------

    m.init();

    // override the region ids
    for ( auto c : m.cells() )
      c->region() = region_ids[c.id()];
    m.set_num_regions( num_elem_blk );

    // creating fields from an exodus file will be problematic since the type
    // information on the file has been thrown away. For example, an int field
    // is written to the file as a floating point field, and there is no way to
    // recover the fact that it is an int field. Furthermore, vector fields
    // are stored as individual scalar fields and recovering the fact that a
    // a collection of scalar fields makes a vector field would require clunky
    // processing of the variable names.


    return status;

#else

    std::cerr << "FLECSI not build with exodus support." << std::endl;
    std::exit(1);

    return -1;

#endif
  }

  //============================================================================
  //!  Implementation of exodus mesh write for burton specialization.
  //!
  //!  \param[in] name Write burton mesh \e m to \e name.
  //!  \param[in] m Burton mesh to write to \e name.
  //!
  //!  \return Exodus error code. 0 on success.
  //============================================================================

  //FIXME: should allow for const mesh_t &
  //int32_t io_exodus_t::write(
  //    const std::string &name, const mesh_t &m) {
  int32_t write( const std::string &name, burton_mesh_t &m ) 
  {

#ifdef HAVE_EXODUS

    std::cout << "Writing mesh to: " << name << std::endl;

    // alias some types
    using std::vector;
    using std::array;
    using   mesh_t = burton_mesh_t;
    using   real_t = typename mesh_t::real_t;
    using integer_t= typename mesh_t::integer_t;
    using vector_t = typename mesh_t::vector_t;
    using ex_real_t = real_t;
    using ex_index_t = int;


    //--------------------------------------------------------------------------
    // initial setup
    //--------------------------------------------------------------------------

    // size of floating point variables used in app.
    int CPU_word_size = sizeof(real_t);
    // size of floating point to be stored in file.
    int IO_word_size = sizeof(ex_real_t);
    // determine the file creation mode
    int cmode = EX_CLOBBER;

    auto exoid =
      ex_create(name.c_str(), cmode, &CPU_word_size, &IO_word_size);
    assert(exoid >= 0);

    // get the general statistics
    auto num_dims = m.num_dimensions();
    auto num_nodes = m.num_vertices();
    auto num_elem = m.num_cells();
    auto num_elem_blk = m.num_regions();
    auto num_node_sets = 0;
    auto num_side_sets = 0;

    // initialize the file.
    auto status = ex_put_init(exoid, "Exodus II output from flecsi.", num_dims,
                              num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets);
    assert(status == 0);

    // get the coordinates from the mesh.
    vector<ex_real_t> xcoord(num_nodes);
    vector<ex_real_t> ycoord(num_nodes);
    for (auto v : m.vertices()) {
      xcoord[v.id()] = v->coordinates()[0];
      ycoord[v.id()] = v->coordinates()[1];
    } // for
    // write the coordinates to the file
    status = ex_put_coord(exoid, xcoord.data(), ycoord.data(), nullptr);
    assert(status == 0);

    // write the coordinate names
    const char *coord_names[3];
    coord_names[0] = "x";
    coord_names[1] = "y";
    coord_names[2] = "z";
    status = ex_put_coord_names(exoid, (char **)coord_names);
    assert(status == 0);





    //--------------------------------------------------------------------------
    // Block connectivity
    //--------------------------------------------------------------------------

    // get the regions
    auto region_cells = m.regions();


    // loop over element blocks
    for ( auto iblk=0; iblk<num_elem_blk; iblk++ ) {

      // set the block header
      auto elem_blk_id = iblk+1;

      // get the elements in this block
      const auto & elem_this_blk = region_cells[iblk];
      auto num_elem_this_blk = elem_this_blk.size();
            
      // count how many vertices there are in this block
      auto num_nodes_this_blk = 0;
      for ( auto c : elem_this_blk ) {
        auto verts = m.vertices(c);
        num_nodes_this_blk += verts.size();
      }

      // set the block header
      auto num_attr_per_elem = 0;
      auto num_faces_per_elem = 0;
      auto num_edges_per_elem = 0;
      status = ex_put_block( 
        exoid, EX_ELEM_BLOCK, elem_blk_id, "nsided", num_elem_this_blk, 
        num_nodes_this_blk, num_edges_per_elem, num_faces_per_elem, 
        num_attr_per_elem
      );
      assert(status == 0);

      // build the connectivitiy list
      vector<ex_index_t> elem_nodes( num_nodes_this_blk );
      vector<ex_index_t> elem_node_counts( num_elem_this_blk );

      auto f = 0, i = 0;
      for ( auto c : elem_this_blk ) {
        auto verts = m.vertices(c);
        elem_node_counts[f++] = verts.size();
        for (auto v : verts)
          elem_nodes[i++] = v.id() + 1; // 1-based ids      
      }
       

      // write connectivity
      status = ex_put_conn(
        exoid, EX_ELEM_BLOCK, elem_blk_id, elem_nodes.data(), 
        nullptr, nullptr
      );
      assert(status == 0);
        
      // write counts
      status = ex_put_entity_count_per_polyhedra(
        exoid, EX_ELEM_BLOCK, elem_blk_id, elem_node_counts.data() 
      );

    } // block

    //--------------------------------------------------------------------------
    // write field data
    //--------------------------------------------------------------------------

    // get the iteration number
    // - first step starts at 1
    auto time_step = 1;


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    //--------------------------------------------------------------------------
    // nodal field data header

    int num_nf = 0; // number of nodal fields
    // real scalars persistent at vertices
    auto rspav = access_type_if(m, real_t, is_persistent_at(vertices));
    num_nf += rspav.size();
    // int scalars persistent at vertices
    auto ispav = access_type_if(m, integer_t, is_persistent_at(vertices));
    num_nf += ispav.size();
    // real vectors persistent at vertices
    auto rvpav = access_type_if(m, vector_t, is_persistent_at(vertices));
    num_nf += num_dims*rvpav.size();

    // variable extension for vectors
    std::string var_ext[3];
    var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

    // put the number of nodal fields
    if (num_nf > 0) {
      status = ex_put_var_param(exoid, "n", num_nf);
      assert(status == 0);
    }

    // fill node variable names array
    int inum = 1;
    for(auto sf: rspav) {
      auto label = validate_string( sf.label() );      
      status = ex_put_var_name(exoid, "n", inum++, label.c_str());
      assert(status == 0);
    } // for
    for(auto sf: ispav) {
      auto label = validate_string( sf.label() );
      status = ex_put_var_name(exoid, "n", inum++, label.c_str());
    } // for
    for(auto vf: rvpav) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        status = ex_put_var_name(exoid, "n", inum++, dim_label.c_str());
      } // for
    } // for


    //--------------------------------------------------------------------------
    // nodal field data

    // write nodal fields
    inum = 1;
    // node field buffer
    vector<ex_real_t> nf(num_nodes);
    for(auto sf: rspav) {
      for(auto v: m.vertices()) nf[v.id()] = sf[v];
      status = ex_put_nodal_var(exoid, time_step, inum++, num_nodes, nf.data());
      assert(status == 0);
    } // for
    for(auto sf: ispav) {
      // cast int fields to real_t
      for(auto v: m.vertices()) nf[v.id()] = (real_t)sf[v];
      status = ex_put_nodal_var(exoid, time_step, inum++, num_nodes, nf.data());
      assert(status == 0);
    } // for
    for(auto vf: rvpav) {
      for(int d=0; d < num_dims; ++d) {
        for(auto v: m.vertices()) nf[v.id()] = vf[v][d];
        status = ex_put_nodal_var(exoid, time_step, inum++, num_nodes, nf.data());
        assert(status == 0);
      } // for
    } // for



    //--------------------------------------------------------------------------
    // element field data headers

    int num_ef = 0; // number of element fields
    // real scalars persistent at cells
    auto rspac = access_type_if(m, real_t, is_persistent_at(cells));
    num_ef += rspac.size();
    // int scalars persistent at cells
    auto ispac = access_type_if(m, integer_t, is_persistent_at(cells));
    num_ef += ispac.size();
    // real vectors persistent at cells
    auto rvpac = access_type_if(m, vector_t, is_persistent_at(cells));
    num_ef += num_dims*rvpac.size();


    // put the number of element fields
    if(num_ef > 0) {
      status = ex_put_var_param(exoid, "e", num_ef);
      assert(status == 0);
    } // if

    // fill element variable names array
    inum = 1;
    for(auto sf: rspac) {
      auto label = validate_string( sf.label() );
      status = ex_put_var_name(exoid, "e", inum++, label.c_str());
      assert(status == 0);
    } // for
    for(auto sf: ispac) {
      auto label = validate_string( sf.label() );
      status = ex_put_var_name(exoid, "e", inum++, label.c_str());
      assert(status == 0);
    } // for
    for(auto vf: rvpac) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        status = ex_put_var_name(exoid, "e", inum++, dim_label.c_str());
      } // for
    } // for


    //--------------------------------------------------------------------------
    // element field data

    // loop over element blocks, writing the filtered element data
    for ( auto iblk=0; iblk<num_elem_blk; iblk++ ) {

      // stats for this block
      auto elem_blk_id = iblk+1;

      // get the elements in this block
      const auto & elem_this_blk = region_cells[iblk];
      auto num_elem_this_blk = elem_this_blk.size();
  
      // write element fields
      inum = 1;
      // element field buffer
      vector<ex_real_t> ef(num_elem_this_blk);
      for(auto sf: rspac) {
        auto cid = 0;
        for(auto c: elem_this_blk) ef[cid++] = sf[c];
        status = ex_put_elem_var(exoid, time_step, inum++, elem_blk_id, num_elem_this_blk, ef.data());
        assert(status == 0);
      } // for
      for(auto sf: ispac) {
        // cast int fields to real_t
        auto cid = 0;
        for(auto c: elem_this_blk) ef[cid++] = (real_t)sf[c];
        status = ex_put_elem_var(exoid, time_step, inum++, elem_blk_id, num_elem_this_blk, ef.data());
        assert(status == 0);
      } // for
      for(auto vf: rvpac) {
        for(int d=0; d < num_dims; ++d) {
          auto cid = 0;
          for(auto c: elem_this_blk) ef[cid++] = vf[c][d];
          status = ex_put_elem_var(exoid, time_step, inum++, elem_blk_id, num_elem_this_blk, ef.data());
          assert(status == 0);
        } // for
      } // for

    } // block

    //--------------------------------------------------------------------------
    // final setup
    //--------------------------------------------------------------------------

    // set the time
    ex_real_t soln_time = m.time();
    status = ex_put_time(exoid, time_step, &soln_time );
    assert(status == 0);

    // close
    status = ex_close(exoid);

    return status;


#else

    std::cerr << "FLECSI not build with exodus support." << std::endl;
    std::exit(1);

    return -1;

#endif

  } // io_exodus_t::write


}; // struct io_exodus_t

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_exodus_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_exodus_t.
//!
//! \return Pointer to io_base_t base class of io_exodus_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_t> * create_io_exodus()
{
  return new burton_io_exodus_t;
} // create_io_exodus






////////////////////////////////////////////////////////////////////////////////
//! Register file extension "g" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_exodus_g_registered =                 
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType(
    "g", create_io_exodus );

////////////////////////////////////////////////////////////////////////////////
//! Register file extension "exo" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_exodus_exo_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType(
    "exo", create_io_exodus );

} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
