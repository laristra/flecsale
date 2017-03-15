/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file 
/// \brief Defines the functionality for the exodus writer and reader.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsi/io/io_base.h"
#include "mesh/burton/burton_mesh.h"


#ifdef HAVE_EXODUS
#  include <exodusII.h>
#endif

// system includes
#include <cstring>

// Paraview has a problem with regions in nfaced data.  Uncomment the next
// line, or compile with -dPARAVIEW_EXODUS_3D_REGION_BUGFIX to outout exodus
// files with only one region.
// #define PARAVIEW_EXODUS_3D_REGION_BUGFIX


namespace ale {
namespace mesh {


////////////////////////////////////////////////////////////////////////////////
/// \brief provides base functionality for exodus writer
/// \tparam N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_io_exodus_base {

public:

  //============================================================================
  //! Typedefs
  //============================================================================

  //! the mesh type
  using mesh_t = burton_mesh_t<N>;

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


#ifdef HAVE_EXODUS

  //============================================================================
  //! \brief open the file for reading or writing
  //! \param [in] name  The name of the file to open.
  //! \param [in] mode  The mode to open the file in.
  //! \return The exodus handle for the open file.
  //============================================================================
  auto open( const std::string &name, std::ios_base::openmode mode ) 
  {
   

#ifdef DEBUG
      // useful for debug
      ex_opts (EX_ABORT | EX_VERBOSE);
#endif

    // size of floating point variables used in app.
    int CPU_word_size = sizeof(real_t);

    if ( (mode & std::ios_base::in) == std::ios_base::in ) 
    {

      // size of floating point stored in name.
      int IO_word_size = 0;
      // the version number
      float version;
      
      // open the file
      exoid_ = ex_open(
        name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version);
      assert(exoid_ >= 0);
      
    }
    else if ( (mode & std::ios_base::out) == std::ios_base::out ) 
    {
      
      // size of floating point to be stored in file.
      int IO_word_size = sizeof(ex_real_t);

      // determine the file creation mode
      int cmode = EX_CLOBBER;

      // create file
      exoid_ =
        ex_create(name.c_str(), cmode, &CPU_word_size, &IO_word_size);
      assert(exoid_ >= 0);

    }

    return exoid_;
  }


  //============================================================================
  //! \brief close the file once completed reading or writing
  //============================================================================
  auto close() 
  {
    auto status = ex_close(exoid_);
    return status;   
  }


  //============================================================================
  //! \brief write the coordinates of the mesh to file.
  //! \return the status of the file
  //============================================================================
  auto write_point_coords( mesh_t & m ) 
  { 
    
    // mesh statistics
    constexpr auto num_dims  = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();

    // storage for coordinates
    std::vector<ex_real_t> coord(num_nodes * num_dims);    

    // copy the coordinates
    for (auto v : m.vertices()) {
      auto & coords = v->coordinates();
      for ( int i=0; i<num_dims; i++ ) 
        coord[ i*num_nodes + v.id() ] = coords[i];
    } // for

    // write the coordinates to the file
    auto status = ex_put_coord( 
      exoid_, coord.data(), coord.data()+num_nodes, coord.data()+2*num_nodes );
    assert(status == 0);
  
    // write the coordinate names
    const char *coord_names[3];
    coord_names[0] = "x";
    coord_names[1] = "y";
    coord_names[2] = "z";
    status = ex_put_coord_names(exoid_, (char **)coord_names);
    assert(status == 0);

    return status;

  }

  //============================================================================
  //! \brief read the coordinates of the mesh from a file.
  //! \return the status of the file
  //============================================================================
  auto read_point_coords( 
    mesh_t & m, size_t num_nodes, std::vector<vertex_t *> & vs ) 
  { 
    // mesh statistics
    constexpr auto num_dims = mesh_t::num_dimensions;

    // intitialize the number of nodes
    m.init_parameters(num_nodes);
    
    // storage for nodes
    vs.clear();
    vs.reserve( num_nodes );

    // read nodes
    std::vector<ex_real_t> coord( num_dims * num_nodes );
    
    auto status = ex_get_coord( 
      exoid_, coord.data(), coord.data()+num_nodes, coord.data()+2*num_nodes);
    assert(status == 0);

    // put nodes into mesh
    for (counter_t i = 0; i < num_nodes; ++i) {
      // convert the point
      point_t p;
      for ( int d=0; d<num_dims; d++ ) 
        p[d] = static_cast<real_t>( coord[ d*num_nodes + i ] );
      // now create it
      auto v = m.create_vertex( p );
      vs.push_back(v);
    } // for

    return status;

  }
    


  //============================================================================
  //! \brief write field data to the file
  //! \param [in] m  The mesh to extract field data from.
  //! \return the status of the file
  //============================================================================
  auto write_fields( mesh_t & m ) 
  { 

    int status;

    // mesh statistics
    constexpr auto num_dims  = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();
    auto num_elem = m.num_cells();
    auto num_elem_blk = m.num_regions();

    //--------------------------------------------------------------------------
    // initial setup

    // get the iteration number
    // - first step starts at 1
    auto time_step = 1;


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

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

    return status;

  }

  //============================================================================
  //! Private Data
  //============================================================================
private:


  //!> the file id
  int exoid_ = -1;

#endif

};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the general templated type for writing and reading from 
///        exodus files.
///
/// \tparam N  The number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_io_exodus_t {};


////////////////////////////////////////////////////////////////////////////////
/// \brief This is the two-dimensional mesh reader and writer based on the 
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_io_exodus_t<2> : 
    public flecsi::io::io_base_t< burton_mesh_t<2> >, burton_io_exodus_base<2>
{

public:

  //============================================================================
  //! \brief Default constructor
  //============================================================================
  burton_io_exodus_t() = default;

  //============================================================================
  //! \brief Implementation of exodus mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  int read( const std::string &name, mesh_t &m) override
  {

#ifdef HAVE_EXODUS

    std::cout << "Reading mesh from: " << name << std::endl;

    // alias some types
    using std::vector;

    // open the exodus file
    auto exoid = open( name, std::ios_base::in );
    assert(exoid >= 0);

    // get the initialization parameters
    ex_init_params exopar;
    auto status = ex_get_init_ext(exoid, &exopar);
    assert(status == 0);

    // verify mesh dimension
    assert(m.num_dimensions == exopar.num_dim);
    auto num_dims = exopar.num_dim;
    auto num_nodes = exopar.num_nodes;
    auto num_elem = exopar.num_elem;
    auto num_elem_blk = exopar.num_elem_blk;
    auto num_node_sets = exopar.num_node_sets;
    auto num_side_sets = exopar.num_side_sets;

    //--------------------------------------------------------------------------
    // read coordinates

    std::vector<vertex_t *> vs;
    status = read_point_coords( m, num_nodes, vs );
    assert( status == 0 );
    

    //--------------------------------------------------------------------------
    // read blocks
    vector<ex_index_t> blockids( num_elem_blk );
    status = ex_get_elem_blk_ids(exoid, blockids.data() );
    assert(status == 0);

    // storage for regions
    vector<ex_index_t> region_ids;
    region_ids.reserve( num_elem );

    // read each block
    for ( int iblk=0; iblk<exopar.num_elem_blk; iblk++ ) {

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
        size_t base = 0;
        for (counter_t e = 0; e < num_elem_this_blk; ++e) {
          elem_vs.clear();
          // get the number of nodes
          num_nodes_per_elem = elem_node_counts[e];
          // copy local vertices into vector ( exodus uses 1 indexed arrays )
          for ( int v=0;  v<num_nodes_per_elem; v++ )
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
        for (counter_t e = 0; e < num_elem_this_blk; ++e) {
          elem_vs.clear();
          // base offset into elt_conn
          auto b = e*num_nodes_per_elem;
          // copy local vertices into vector ( exodus uses 1 indexed arrays )
          for ( int v=0;  v<num_nodes_per_elem; v++ )
            elem_vs.emplace_back( vs[ elt_conn[b+v]-1 ] );
          // create acual cell
          auto c = m.create_cell( elem_vs );          
        }

      } // element type
      //--------------------------------


      // set element regions
      for ( counter_t e = 0; e < num_elem_this_blk; e++ )
        region_ids.emplace_back( iblk );
      
    }
    // end blocks
    //--------------------------------------------------------------------------

    m.init();

    // override the region ids
    m.set_regions( region_ids.data() );
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
  //!  \brief Implementation of exodus mesh write for burton specialization.
  //!
  //!  \param[in] name Write burton mesh \e m to \e name.
  //!  \param[in] m Burton mesh to write to \e name.
  //!
  //!  \return Exodus error code. 0 on success.
  //============================================================================

  //FIXME: should allow for const mesh_t &
  //int io_exodus_t::write(
  //    const std::string &name, const mesh_t &m) {
  int write( const std::string &name, mesh_t &m ) override
  {

#ifdef HAVE_EXODUS

    std::cout << "Writing mesh to: " << name << std::endl;

    // alias some types
    using std::vector;
    using std::array;


    //--------------------------------------------------------------------------
    // initial setup
    //--------------------------------------------------------------------------

    auto exoid = open( name, std::ios_base::out );
    assert(exoid >= 0);

    // get the general statistics
    constexpr auto num_dims = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();
    auto num_elem = m.num_cells();
    auto num_elem_blk = m.num_regions();
    auto num_node_sets = 0;
    auto num_side_sets = 0;

    // initialize the file.
    auto status = ex_put_init(exoid, "Exodus II output from flecsi.", num_dims,
                              num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets);
    assert(status == 0);

    //--------------------------------------------------------------------------
    // Point Coordinates
    //--------------------------------------------------------------------------
    
    status = write_point_coords( m );
    assert( status == 0 );


    //--------------------------------------------------------------------------
    // Block connectivity
    //--------------------------------------------------------------------------

    // get the regions
    auto region_cells = m.regions();


    // loop over element blocks
    for ( int iblk=0; iblk<num_elem_blk; iblk++ ) {

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

      size_t f = 0, i = 0;
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
    status = write_fields( m );
    assert( status == 0 );


    //--------------------------------------------------------------------------
    // final step
    //--------------------------------------------------------------------------
    status = close();

    return status;


#else

    std::cerr << "FLECSI not build with exodus support." << std::endl;

    return 0;

#endif

  } // io_exodus_t::write


}; // struct io_exodus_t






////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the 
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_io_exodus_t<3> : 
    public flecsi::io::io_base_t< burton_mesh_t<3> >, burton_io_exodus_base<3>

{

public:

  //============================================================================
  //! \brief Default constructor
  //============================================================================
  burton_io_exodus_t() = default;


  //============================================================================
  //! \brief Implementation of exodus mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  int read( const std::string &name, mesh_t &m) override
  {

#ifdef HAVE_EXODUS

    std::cout << "Reading mesh from: " << name << std::endl;

    // alias some types
    using std::vector;

    // open the exodus file
    auto exoid = open( name, std::ios_base::in );

    // get the initialization parameters
    ex_init_params exopar;
    auto status = ex_get_init_ext(exoid, &exopar);
    assert(status == 0);

    // verify mesh dimension
    assert(m.num_dimensions == exopar.num_dim);
    auto num_dims = exopar.num_dim;
    auto num_nodes = exopar.num_nodes;
    auto num_elem = exopar.num_elem;
    auto num_elem_blk = exopar.num_elem_blk;
    auto num_face_blk = exopar.num_face_blk;
    auto num_node_sets = exopar.num_node_sets;
    auto num_side_sets = exopar.num_side_sets;

    //--------------------------------------------------------------------------
    // read coordinates

    std::vector<vertex_t *> vertices;
    status = read_point_coords( m, num_nodes, vertices );
    assert( status == 0 );

    //--------------------------------------------------------------------------
    // read face blocks

    vector<ex_index_t> face_block_ids( num_face_blk );

    if ( num_face_blk > 0 ) {
      status = ex_get_ids(exoid, EX_FACE_BLOCK, face_block_ids.data() );
      assert(status == 0);
    }

    // storage for faces
    std::vector<face_t *> faces;

    // read each block
    for ( int iblk=0; iblk<num_face_blk; iblk++ ) {

      auto face_blk_id = face_block_ids[iblk];

      // get the info about this block
      ex_index_t num_face_this_blk = 0;
      ex_index_t num_faces_per_face = 0;
      ex_index_t num_edges_per_face = 0;
      ex_index_t num_nodes_per_face = 0;
      ex_index_t num_attr = 0;
      char face_type[MAX_STR_LENGTH];

      status = ex_get_block(
        exoid, EX_FACE_BLOCK, face_blk_id, face_type, &num_face_this_blk, 
        &num_nodes_per_face, &num_edges_per_face, &num_faces_per_face, &num_attr );
      assert(status == 0);
      

      // resize face storage for new faces
      faces.reserve( faces.size() + num_face_this_blk );

      // the number of nodes per face is really the number of
      // nodes in the whole block ( includes duplicate / overlapping
      // nodes )
      auto num_nodes_this_blk = num_nodes_per_face;

      // get the number of nodes per element
      vector<ex_index_t> face_node_counts( num_nodes_this_blk );
      status = ex_get_entity_count_per_polyhedra( 
        exoid, EX_FACE_BLOCK, face_blk_id, face_node_counts.data() );
      assert(status == 0);
      
      // read element definitions
      vector<ex_index_t> face_nodes( num_nodes_this_blk );
      status = ex_get_conn( 
        exoid, EX_FACE_BLOCK, face_blk_id, face_nodes.data(), nullptr, nullptr );
      assert(status == 0);
        
      // storage for element verts
      vector<vertex_t *> face_vs;
      face_vs.reserve( face_node_counts[0] );
        
      // create faces in mesh
      for (counter_t e=0, base=0; e < num_face_this_blk; ++e) {
        face_vs.clear();
        // get the number of nodes
        num_nodes_per_face = face_node_counts[e];
        // copy local vertices into vector ( exodus uses 1 indexed arrays )
        for ( auto v=0;  v<num_nodes_per_face; v++ )
          face_vs.emplace_back( vertices[ face_nodes[base+v] - 1 ] );
        // create acual face
        auto f = m.create_face( face_vs );
        // now add it to the 
        faces.emplace_back( std::move(f) );
        // base offset into face conn
        base += num_nodes_per_face;
      }

    }

    //--------------------------------------------------------------------------
    // read element blocks
    vector<ex_index_t> elem_block_ids( num_elem_blk );
    status = ex_get_elem_blk_ids(exoid, elem_block_ids.data() );
    assert(status == 0);

    // storage for regions
    vector<ex_index_t> region_ids;
    region_ids.reserve( num_elem );

    // need to keep track of who owns each face
    // this is needed to ascertain a face direction
    vector< std::pair< face_t*, cell_t* > > 
      face_owner( faces.size(), std::make_pair(nullptr, nullptr)  );

    // read each block
    for ( int iblk=0; iblk<num_elem_blk; iblk++ ) {

      auto elem_blk_id = elem_block_ids[iblk];

      char block_name[MAX_LINE_LENGTH];
      status = ex_get_name(exoid, EX_ELEM_BLOCK, elem_blk_id, block_name);
      assert(status == 0);

      // get the info about this block
      ex_index_t num_elem_this_blk = 0;
      ex_index_t num_faces_per_elem = 0;
      ex_index_t num_edges_per_elem = 0;
      ex_index_t num_nodes_per_elem = 0;
      ex_index_t num_attr = 0;
      char elem_type[MAX_STR_LENGTH];
      status = ex_get_block(
        exoid, EX_ELEM_BLOCK, elem_blk_id, elem_type, &num_elem_this_blk, 
        &num_nodes_per_elem, &num_edges_per_elem, &num_faces_per_elem, &num_attr );
      assert(status == 0);

      //--------------------------------
      // polygon data
      if ( strcmp("nfaced",elem_type) == 0 || strcmp("NFACED",elem_type) == 0 ) {

        // the number of faces per element is really the number of
        // faces in the whole block ( includes duplicate / overlapping
        // nodes )
        auto num_face_this_blk = num_faces_per_elem;

        // get the number of nodes per element
        vector<ex_index_t> elem_face_counts( num_face_this_blk );
        status = ex_get_entity_count_per_polyhedra( 
          exoid, EX_ELEM_BLOCK, elem_blk_id, elem_face_counts.data() );
        assert(status == 0);

        // read element definitions
        vector<ex_index_t> elem_faces( num_face_this_blk );
        status = ex_get_conn( 
          exoid, EX_ELEM_BLOCK, elem_blk_id, nullptr, nullptr, elem_faces.data() );
        assert(status == 0);
        
        // storage for element faces
        vector<face_t *>   elem_fs;
        vector<ex_index_t> elem_fs_ids;
        elem_fs.reserve( elem_face_counts[0] );
        elem_fs_ids.reserve( elem_face_counts[0] );
        
        // create cells in mesh
        for (counter_t e=0, base=0; e < num_elem_this_blk; ++e) {
          // reset storage
          elem_fs.clear();
          // get the number of faces
          num_faces_per_elem = elem_face_counts[e];
          // copy local vertices into vector ( exodus uses 1 indexed arrays )
          for ( int v=0;  v<num_faces_per_elem; v++ ) {
            auto id = elem_faces[base+v] - 1;            
            elem_fs.emplace_back( faces[ id ] );
            elem_fs_ids.emplace_back( id );
          }
          // create acual cell
          auto c = m.create_cell( elem_fs );
          // Since exodus does not provide the face directions, we need to figure 
          // these out with geometric checks
          for ( auto face_id : elem_fs_ids ) {
            if ( ! face_owner[face_id].first ) {
              auto f = faces[face_id];
              face_owner[face_id].first  = f;
              face_owner[face_id].second = c;
            }
          }
          // base offset into elt_conn
          base += num_faces_per_elem;
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
        for (counter_t e = 0; e < num_elem_this_blk; ++e) {
          elem_vs.clear();
          // base offset into elt_conn
          auto b = e*num_nodes_per_elem;
          // copy local vertices into vector ( exodus uses 1 indexed arrays )
          for ( int v=0;  v<num_nodes_per_elem; v++ ) 
            elem_vs.emplace_back( vertices[ elt_conn[b+v]-1 ] );
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
    // final mesh setup
    m.init();

#ifndef PARAVIEW_EXODUS_3D_REGION_BUGFIX

    // override the region ids
    m.set_regions( region_ids.data() );
    m.set_num_regions( num_elem_blk );

#endif

    // loop over faces that got created in the face blocks and and
    // make sure the owner is the first cell
    for ( auto face_pair : face_owner ) {
      auto face = face_pair.first;
      auto cell = face_pair.second;
      auto cells = m.cells( face );
      assert( cells.front().id() == cell->template id<0>() );
    }
    


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
  //!  \brief Implementation of exodus mesh write for burton specialization.
  //!
  //!  \param[in] name Write burton mesh \e m to \e name.
  //!  \param[in] m Burton mesh to write to \e name.
  //!
  //!  \return Exodus error code. 0 on success.
  //============================================================================

  //FIXME: should allow for const mesh_t &
  //int io_exodus_t::write(
  //    const std::string &name, const mesah_t &m) {
  int write( const std::string &name, mesh_t &m ) override
  {

#ifdef HAVE_EXODUS

    std::cout << "Writing mesh to: " << name << std::endl;

    // alias some types
    using std::vector;
    using std::array;


    //--------------------------------------------------------------------------
    // initial setup
    //--------------------------------------------------------------------------

    auto exoid = open( name, std::ios_base::out );
    assert(exoid >= 0);


    // get the general statistics
    constexpr auto num_dims = mesh_t::num_dimensions;
    auto num_nodes = m.num_vertices();
    auto num_faces = m.num_faces();
    auto num_elem = m.num_cells();
    auto num_elem_blk = m.num_regions();

    // set exodus parameters
    ex_init_params exopar;
    strcpy( exopar.title, "Exodus II output from flecsi." );
    exopar.num_dim = num_dims;
    exopar.num_nodes = num_nodes;
    exopar.num_edge = 0;
    exopar.num_edge_blk = 0;
    exopar.num_face = num_faces;
    exopar.num_face_blk = 1;
    exopar.num_elem = num_elem;
    exopar.num_elem_blk = num_elem_blk;
    exopar.num_node_sets = 0;
    exopar.num_edge_sets = 0;
    exopar.num_face_sets = 0;
    exopar.num_side_sets = 0;
    exopar.num_elem_sets = 0;    
    exopar.num_node_maps = 0;
    exopar.num_edge_maps = 0;
    exopar.num_face_maps = 0;
    exopar.num_elem_maps = 0;
    auto status = ex_put_init_ext( exoid, &exopar );
    assert(status == 0);

    //--------------------------------------------------------------------------
    // Point Coordinates
    //--------------------------------------------------------------------------
    
    status = write_point_coords( m );
    assert( status == 0 );



    //--------------------------------------------------------------------------
    // Face Connectivity
    //--------------------------------------------------------------------------
    {

      // set the block header
      auto face_blk_id = num_elem_blk + 1; // don't collide
      
      // put all faces in one block
      auto faces_this_blk = m.faces();
      auto num_faces_this_blk = faces_this_blk.size();

      // count how many vertices and faces there are in this block
      // double count face nodes
      auto num_nodes_this_blk = 0;
      for ( auto f : faces_this_blk )
        num_nodes_this_blk += m.vertices(f).size();
      
      // set the block header
      auto num_attr_per_face = 0;
      auto num_faces_per_face = 0;
      auto num_edges_per_face = 0;
      status = ex_put_block( 
        exoid, EX_FACE_BLOCK, face_blk_id, "nsided", num_faces_this_blk, 
        num_nodes_this_blk, num_edges_per_face, num_faces_per_face, 
        num_attr_per_face
      );
      assert(status == 0);
      
      // the block name
      char name[256];
      sprintf( name, "face_block_%lu", face_blk_id );
      status = ex_put_name( exoid, EX_FACE_BLOCK, face_blk_id, name );

      // build the connectivitiy list
      vector<ex_index_t> face_nodes( num_nodes_this_blk );
      vector<ex_index_t> face_node_counts( num_faces_this_blk );

      size_t e = 0, i = 0;
      // for each face, get nodes
      for ( auto f : faces_this_blk ) {
        // node count
        auto verts = m.vertices(f);
        face_node_counts[e++] = verts.size();
        // vertex ids
        for (auto v : verts )
          face_nodes[i++] = v.id() + 1; // 1-based ids      
      }
       
      // write connectivity
      status = ex_put_conn(
        exoid, EX_FACE_BLOCK, face_blk_id, face_nodes.data(), 
        nullptr, nullptr
      );
      assert(status == 0);
        
      // write counts
      status = ex_put_entity_count_per_polyhedra(
        exoid, EX_FACE_BLOCK, face_blk_id, face_node_counts.data() 
      );
      assert(status == 0);

    } // face block

    //--------------------------------------------------------------------------
    // Element Block connectivity
    //--------------------------------------------------------------------------

    // get the regions
    auto region_cells = m.regions();

    // loop over element blocks
    for ( int iblk=0; iblk<num_elem_blk; iblk++ ) {

      // set the block header
      auto elem_blk_id = iblk+1;

      // get the elements in this block
      const auto & elem_this_blk = region_cells[iblk];
      auto num_elem_this_blk = elem_this_blk.size();

      // count how many faces there are in this block
      // ( double count faces )
      auto num_faces_this_blk = 0;
      for ( auto c : elem_this_blk )
        num_faces_this_blk += m.faces(c).size();

      // set the block header
      auto num_attr_per_elem = 0;
      auto num_edges_per_elem = 0;
      auto num_nodes_per_elem = 0;
      status = ex_put_block( 
        exoid, EX_ELEM_BLOCK, elem_blk_id, "nfaced", num_elem_this_blk, 
        num_nodes_per_elem, num_edges_per_elem, num_faces_this_blk, 
        num_attr_per_elem
      );
      assert(status == 0);

      // the block name
      char name[256];
      sprintf( name, "elem_block_%d", elem_blk_id );
      status = ex_put_name( exoid, EX_ELEM_BLOCK, elem_blk_id, name );
      assert(status == 0);

      // build the connectivitiy list
      vector<ex_index_t> elem_faces( num_faces_this_blk );
      vector<ex_index_t> elem_face_counts( num_elem_this_blk );

      size_t e = 0, i = 0;
      for ( auto c : elem_this_blk ) {
        // get faces of this element
        auto faces = m.faces(c);
        // face count
        elem_face_counts[e++] = faces.size();
        // for each element face, get ids
        for ( auto f : faces ) 
          elem_faces[i++] = f.id() + 1; // 1-based ids      
      }
       
      // write connectivity
      status = ex_put_conn(
        exoid, EX_ELEM_BLOCK, elem_blk_id, nullptr, nullptr,
        elem_faces.data() 
      );
      assert(status == 0);
        
      // write counts
      status = ex_put_entity_count_per_polyhedra(
        exoid, EX_ELEM_BLOCK, elem_blk_id, elem_face_counts.data() 
      );
      assert(status == 0);
            
    } // block

    //--------------------------------------------------------------------------
    // write field data
    //--------------------------------------------------------------------------
    status = write_fields( m );
    assert( status == 0 );


    //--------------------------------------------------------------------------
    // final step
    //--------------------------------------------------------------------------
    status = close();

    return status;


#else

    std::cerr << "FLECSI not build with exodus support." << std::endl;

    return 0;

#endif

  } // io_exodus_t::write



private:

  //! \brief select an appropriate mapping function
  static 
  std::function<unsigned int (unsigned int)> get_face_mapper( size_t num_faces ) 
  {
    //! \brief the exodus-to-internal face mapping for hexes
    constexpr unsigned int hex_face_map_[] = { 2, 3, 4, 5, 0, 1 };
    //! \brief the exodus-to-internal face mapping for tets
    constexpr unsigned int tet_face_map_[] = { 2, 0, 1, 3 };
    
    switch( num_faces ) {
    case 4:
      return [&]( unsigned int i ) { return tet_face_map_[i]; };
    case 6:
      return [&]( unsigned int i ) { return hex_face_map_[i]; };
    default:
      return []( unsigned int i ) { return i; };
    }
  }


}; // struct io_exodus_t



////////////////////////////////////////////////////////////////////////////////
//! \brief Create a burton_io_exodus_t and return a pointer to the base class.
//!
//! \tparam N  The number of mesh dimensions.
//!
//! \return Pointer to io_base_t base class of io_exodus_t.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
inline flecsi::io::io_base_t< burton_mesh_t<N> > * create_io_exodus()
{
  return new burton_io_exodus_t<N>;
} // create_io_exodus



////////////////////////////////////////////////////////////////////////////////
//! Register file extension "g" with factory.
////////////////////////////////////////////////////////////////////////////////
//! @{
static bool burton_2d_exodus_g_registered =                 
  flecsi::io::io_factory_t< burton_mesh_t<2> >::instance().registerType(
    "g", create_io_exodus<2> );

static bool burton_3d_exodus_g_registered =                 
  flecsi::io::io_factory_t< burton_mesh_t<3> >::instance().registerType(
    "g", create_io_exodus<3> );
//! @}

////////////////////////////////////////////////////////////////////////////////
//! Register file extension "exo" with factory.
////////////////////////////////////////////////////////////////////////////////
//! @{
static bool burton_2d_exodus_exo_registered =
  flecsi::io::io_factory_t< burton_mesh_t<2> >::instance().registerType(
    "exo", create_io_exodus<2> );

static bool burton_3d_exodus_exo_registered =
  flecsi::io::io_factory_t< burton_mesh_t<3> >::instance().registerType(
    "exo", create_io_exodus<3> );
//! @}

} // namespace mesh
} // namespace ale
