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
#include <fstream>


#ifdef HAVE_TECIO
#  include <TECIO.h>
#endif

//! user includes
#include "flecsi/io/io_base.h"
#include "../../mesh/burton/burton_mesh.h"
#include "../../utils/errors.h"
#include "../../utils/string_utils.h"


namespace ale {
namespace mesh {

//! class for variable locations
enum class tec_var_location_t 
{
  cell = 0,
  node = 1
};


////////////////////////////////////////////////////////////////////////////////
/// \class io_tecplot_ascii_t io_tecplot.h
/// \brief io_tecplot_ascii_t provides a derived type of io_base.h and registrations
///   of the tecplot_ascii file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_tecplot_ascii_t : public flecsi::io_base_t<burton_mesh_t> {
  
  //! Default constructor
  burton_io_tecplot_ascii_t() {}

  //============================================================================
  //! Implementation of tecplot mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return tecplot error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_t &m) override
  {


    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    std::cout << "Writing mesh to: " << name << std::endl;

    // some aliazes
    using std::endl;
    using std::pair;
    using std::make_pair;
    using std::string;
    using std::vector;

    using   mesh_t = burton_mesh_t;
    using   real_t = typename mesh_t::real_t;
    using vector_t = typename mesh_t::vector_t;

    // get the general statistics
    auto num_dims  = m.num_dimensions();

  
    // open the file for writing
    std::ofstream ofs( name.c_str() );
    assert( ofs.good() && "error opening file" );

    //--------------------------------------------------------------------------
    // collect field data
    //--------------------------------------------------------------------------


    // variable extension for vectors
    std::string var_ext[3];
    var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

    // collect the variables
    vector< pair<string,tec_var_location_t> > variables;

    // write the coordinate names
    variables.emplace_back( make_pair( "x", tec_var_location_t::node ) );
    variables.emplace_back( make_pair( "y", tec_var_location_t::node ) );


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    //----------------------------------------------------------------------------
    // nodal field data
    //----------------------------------------------------------------------------

    int num_nf = 0; // number of nodal fields
    // real scalars persistent at vertices
    auto rspav = access_type_if(m, real_t, is_persistent_at(vertices));
    num_nf += rspav.size();
    // int scalars persistent at vertices
    auto ispav = access_type_if(m, int, is_persistent_at(vertices));
    num_nf += ispav.size();
    // real vectors persistent at vertices
    auto rvpav = access_type_if(m, vector_t, is_persistent_at(vertices));
    num_nf += num_dims*rvpav.size();

    // fill node variable names array
    for(auto sf: rspav) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::node ) );
    }
    for(auto sf: ispav) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::node ) );
    }
    for(auto vf: rvpav) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        variables.emplace_back( make_pair( dim_label, tec_var_location_t::node ) );
      } // for
    } // for
  

    //----------------------------------------------------------------------------
    // element field data
    //----------------------------------------------------------------------------

    int num_ef = 0; // number of element fields
    // real scalars persistent at cells
    auto rspac = access_type_if(m, real_t, is_persistent_at(cells));
    num_ef += rspac.size();
    // int scalars persistent at cells
    auto ispac = access_type_if(m, int, is_persistent_at(cells));
    num_ef += ispac.size();
    // real vectors persistent at cells
    auto rvpac = access_type_if(m, vector_t, is_persistent_at(cells));
    num_ef += num_dims*rvpac.size();


    // fill element variable names array
    for(auto sf: rspac) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::cell ) );
    }
    for(auto sf: ispac) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::cell ) );
    }
    for(auto vf: rvpac) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        variables.emplace_back( make_pair( dim_label, tec_var_location_t::cell ) );
      } // for
    } // for

    //--------------------------------------------------------------------------
    // HEADER
    //--------------------------------------------------------------------------


    ofs << "TITLE = \"Tecplot output from flecsi.\"" << endl;
    ofs << "FILETYPE = FULL" << endl;

    ofs << "VARIABLES = " << endl;
    for ( const auto & label : variables )
      ofs << "\"" << label.first << "\" ";
    ofs << endl;

    ofs << "ZONE" << endl;
    ofs << "T = \"zone block\" "
        << "ZONETYPE=FEPOLYGON "
        << "NODES=" << m.num_vertices() << " FACES=" << m.num_edges() << " ELEMENTS=" << m.num_cells() << " "
        << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0 DATAPACKING=BLOCK ";

    auto ef_start = num_dims+num_nf+1;
    auto ef_end = num_dims+num_nf+num_ef;
    if ( num_ef > 0 )
      ofs << "VARLOCATION=([" << ef_start << "-" << ef_end << "]=CELLCENTERED)";
    
    ofs << endl;

    //--------------------------------------------------------------------------
    // WRITE DATA
    //--------------------------------------------------------------------------

    // get the coordinates from the mesh.
    for(int d=0; d < num_dims; ++d) {
      for (auto v : m.vertices()) 
        ofs << v->coordinates()[d] << " ";
      ofs << endl;
    } // for

    //----------------------------------------------------------------------------
    // nodal field data
    //----------------------------------------------------------------------------

    // node field buffer
    for(auto sf: rspav) {
      for(auto v: m.vertices()) ofs << sf[v] << " ";
      ofs << endl;
    } // for
    for(auto sf: ispav) {
      for(auto v: m.vertices()) ofs << sf[v] << " ";
      ofs << endl;
    } // for
    for(auto vf: rvpav) {
      for(int d=0; d < num_dims; ++d) {
        for(auto v: m.vertices()) ofs << vf[v][d] << " ";
        ofs << endl;
      } // for
    } // for

    //----------------------------------------------------------------------------
    // nodal field data
    //----------------------------------------------------------------------------

    // element field buffer
    for(auto sf: rspac) {
      for(auto c: m.cells()) ofs << sf[c] << " ";
      ofs << endl;
    } // for
    for(auto sf: ispac) {
      for(auto c: m.cells()) ofs << sf[c] << " ";
      ofs << endl;
    } // for
    for(auto vf: rvpac) {
      for(int d=0; d < num_dims; ++d) {
        for(auto c: m.cells()) ofs << vf[c][d] << " ";
        ofs << endl;
      } // for
    } // for

    //--------------------------------------------------------------------------
    // WRITE CONNECTIVITY
    //--------------------------------------------------------------------------

    // the nodes of each face
    ofs << "#face nodes" << endl;
    for (auto e : m.edges()) {
      for (auto v : m.vertices(e))
        ofs << v.id() + 1 << " "; // 1-based ids
      ofs << endl;
    }

    ofs << "#left elements" << endl;
    for (auto e : m.edges()) {
      auto cells = m.cells( e );
      assert( cells.size() > 0 && "cell has no edges");
      // always has left cell
      ofs << cells[0].id() + 1 << " "; // 1-based ids
    }
    ofs << endl;


    ofs << "#right elements" << endl;
    for (auto e : m.edges()) {
      auto cells = m.cells( e );
      // boundary faces don't have right cell
      auto id = ( cells.size() == 2 ) ? cells[1].id()+1 : 0;
      ofs << id << " ";
    }
    ofs << endl;
 
    //--------------------------------------------------------------------------
    // Finalize
    //--------------------------------------------------------------------------
  
    // close file stream
    ofs.close();

    return 0;

  } // io_tecplot_ascii_t::write


  //============================================================================
  //! Implementation of tecplot mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to Read to \e name.
  //!
  //! \return tecplot error code. 0 on success.
  //!
  //============================================================================
  int32_t read( const std::string &name, burton_mesh_t &m)  override
  {
    raise_implemented_error( "No tecplot read functionality has been implemented" );
  };

}; // io_tecplot_ascii_t

////////////////////////////////////////////////////////////////////////////////
/// \class io_tecplot_binary_t io_tecplot.h
/// \brief io_tecplot_binary_t provides a derived type of io_base.h and registrations
///   of the tecplot_binary file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
struct burton_io_tecplot_binary_t : public flecsi::io_base_t<burton_mesh_t> {

  //! Default constructor
  burton_io_tecplot_binary_t() {}

  //============================================================================
  //! Implementation of tecplot mesh write for burton specialization.
  //!
  //! \param[in] name Write burton mesh \e m to \e name.
  //! \param[in] m Burton mesh to write to \e name.
  //!
  //! \return tecplot error code. 0 on success.
  //!
  //! FIXME: should allow for const mesh_t &
  //============================================================================
  int32_t write( const std::string &name, burton_mesh_t &m)  override
  {

#ifdef HAVE_TECIO

    std::cout << "Writing mesh to: " << name << std::endl;
    

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    // alias some types
    using std::pair;
    using std::make_pair;
    using std::string;
    using std::vector;

    using   mesh_t = burton_mesh_t;
    using   real_t = typename mesh_t::real_t;
    using vector_t = typename mesh_t::vector_t;
    using tec_real_t = real_t;
    using tec_int_t = INTEGER4;

    /* Open the file & write the datafile header information */
    tec_int_t Debug     = 0; // Set to 0 for no debugging or 1 to debug
    tec_int_t FileType  = 0; // 0=Tecplot binary (.plt) 1=Tecplot subzone (.szplt)

    tec_int_t VIsDouble; // 0=Single 1=Double
    if ( std::is_same_v<tec_real_t, float> )
      VIsDouble = 0;
    else if ( std::is_same_v<tec_real_t, double> )
      VIsDouble = 1;
    else
      raise_implemented_error( "Can only output to tecplot with floats or doubls" );

    // get the general statistics
    tec_int_t num_dims  = m.num_dimensions();
    tec_int_t num_nodes = m.num_vertices();
    tec_int_t num_edges = m.num_edges();
    tec_int_t num_faces = num_edges;
    tec_int_t num_elem  = m.num_cells();

    // set the time
    double soln_time = m.get_time();

    // variable extension for vectors
    std::string var_ext[3];
    var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

    // collect the variables
    vector< pair<string,tec_var_location_t> > variables;

    // write the coordinate names
    variables.emplace_back( make_pair( "x", tec_var_location_t::node ) );
    variables.emplace_back( make_pair( "y", tec_var_location_t::node ) );

    //--------------------------------------------------------------------------
    // collect field data
    //--------------------------------------------------------------------------

    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return utils::replace_all( std::forward<decltype(str)>(str), " ", "_" );;
    };



    //----------------------------------------------------------------------------
    // nodal field data

    // real scalars persistent at vertices
    auto rspav = access_type_if(m, real_t, is_persistent_at(vertices));
    // int scalars persistent at vertices
    auto ispav = access_type_if(m, int, is_persistent_at(vertices));
    // real vectors persistent at vertices
    auto rvpav = access_type_if(m, vector_t, is_persistent_at(vertices));

    // fill node variable names array
    for(auto sf: rspav) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::node ) );
    }
    for(auto sf: ispav) {
      auto label = validate_string( sf.label() );      
      variables.emplace_back( make_pair( label, tec_var_location_t::node ) );
    }
    for(auto vf: rvpav) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        variables.emplace_back( make_pair( dim_label, tec_var_location_t::node ) );
      } // for
    } // for
  

    //----------------------------------------------------------------------------
    // element field data

    // real scalars persistent at cells
    auto rspac = access_type_if(m, real_t, is_persistent_at(cells));
    // int scalars persistent at cells
    auto ispac = access_type_if(m, int, is_persistent_at(cells));
    // real vectors persistent at cells
    auto rvpac = access_type_if(m, vector_t, is_persistent_at(cells));


    // fill element variable names array
    for(auto sf: rspac) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::cell ) );
    }
    for(auto sf: ispac) {
      auto label = validate_string( sf.label() );
      variables.emplace_back( make_pair( label, tec_var_location_t::cell ) );
    }
    for(auto vf: rvpac) {
      auto label = validate_string( vf.label() );
      for(int d=0; d < num_dims; ++d) {
        auto dim_label = label + var_ext[d];
        variables.emplace_back( make_pair( dim_label, tec_var_location_t::cell ) );
      } // for
    } // for

    //--------------------------------------------------------------------------
    // WRITE HEADERS
    //--------------------------------------------------------------------------
    
    // create the variable name string
    string var_string;
    for ( const auto & label : variables )
      var_string += label.first + " ";

    auto status = TECINI112( const_cast<char*>( "Tecplot output from flecsi." ),
                             const_cast<char*>( var_string.c_str() ),
                             const_cast<char*>( name.c_str() ),
                             const_cast<char*>( "." ), // scratch dir
                             &FileType,
                             &Debug,
                             &VIsDouble );
    assert( status == 0 && "error with TECINI" );

  
    //----------------------------------------------------------------------------
    // Create ZONE header
    //----------------------------------------------------------------------------


    tec_int_t ZoneType = 6; // set the zone type to FEPolygon
    tec_int_t NumFaces = 0; // not used for for most zone types
    tec_int_t ICellMax = 0; // reserved for future use
    tec_int_t JCellMax = 0; // reserved for future use
    tec_int_t KCellMax = 0; // reserved for future use

    tec_int_t StrandID = 0; // StaticZone
    tec_int_t ParentZn = 0; // no parent
    tec_int_t IsBlock  = 1; // this is a Block
    tec_int_t NumFaceConnections       = 0; // not used 
    tec_int_t FaceNeighborMode         = 0; // not used
    tec_int_t TotalNumFaceNodes        = 2*num_faces; // total nodes for all faces
    tec_int_t TotalNumBndryFaces       = 0; // boundary faces not used
    tec_int_t TotalNumBndryConnections = 0; // boundary face connections not used
    tec_int_t ShareConnectivityFromZone = 0;  // pass 0 to indicate no sharign
    tec_int_t * PassiveVarList = nullptr;
    tec_int_t * ShareVarFromZone = nullptr;

    // create the variable location array
    vector<tec_int_t> var_locations;
    for ( const auto & var : variables )
      var_locations.push_back( static_cast<tec_int_t>(var.second) );

    status = TECZNE112( const_cast<char*>( "zone block" ),
                        &ZoneType,
                        &num_nodes,
                        &num_elem,
                        &num_faces,
                        &ICellMax,
                        &JCellMax,
                        &KCellMax,
                        &soln_time,
                        &StrandID,
                        &ParentZn,
                        &IsBlock,
                        &NumFaceConnections,
                        &FaceNeighborMode,
                        &TotalNumFaceNodes,
                        &TotalNumBndryFaces,
                        &TotalNumBndryConnections,
                        PassiveVarList,
                        var_locations.data(),
                        ShareVarFromZone,
                        &ShareConnectivityFromZone);
  
    assert( status == 0 && "error with TECZNE" );



    //--------------------------------------------------------------------------
    // WRITE DATA
    //--------------------------------------------------------------------------

    // a temporary for vector values
    vector<tec_real_t> xvals( std::max( num_nodes, num_elem ) );
    vector<tec_real_t> yvals( std::max( num_nodes, num_elem ) );

    // get the coordinates from the mesh.
    for (auto v : m.vertices()) {
      xvals[v.id()] = v->coordinates()[0];
      yvals[v.id()] = v->coordinates()[1];
    } // for

      // write the coordinates to the file
    status = TECDAT112( &num_nodes, xvals.data(), &VIsDouble );
    assert( status == 0 && "error with TECDAT" );

    status = TECDAT112( &num_nodes, yvals.data(), &VIsDouble );
    assert( status == 0 && "error with TECDAT" );

    //----------------------------------------------------------------------------
    // nodal field data

    // node field buffer
    for(auto sf: rspav) {
      for(auto v: m.vertices()) xvals[v.id()] = sf[v];
      status = TECDAT112( &num_nodes, xvals.data(), &VIsDouble );
      assert( status == 0 && "error with TECDAT" );
    } // for
    for(auto sf: ispav) {
      // cast int fields to real_t
      for(auto v: m.vertices()) xvals[v.id()] = (tec_real_t)sf[v];
      status = TECDAT112( &num_nodes, xvals.data(), &VIsDouble );
      assert( status == 0 && "error with TECDAT" );
    } // for
    for(auto vf: rvpav) {
      for(int d=0; d < num_dims; ++d) {
        for(auto v: m.vertices()) xvals[v.id()] = vf[v][d];
        status = TECDAT112( &num_nodes, xvals.data(), &VIsDouble );
        assert( status == 0 && "error with TECDAT" );
      } // for
    } // for

    //----------------------------------------------------------------------------
    // cell field data

    // element field buffer
    for(auto sf: rspac) {
      for(auto c: m.cells()) xvals[c.id()] = sf[c];
      status = TECDAT112( &num_elem, xvals.data(), &VIsDouble );
      assert( status == 0 && "error with TECDAT" );
    } // for
    for(auto sf: ispac) {
      // cast int fields to real_t
      for(auto c: m.cells()) xvals[c.id()] = (tec_real_t)sf[c];
      status = TECDAT112( &num_elem, xvals.data(), &VIsDouble );
      assert( status == 0 && "error with TECDAT" );
    } // for
    for(auto vf: rvpac) {
      for(int d=0; d < num_dims; ++d) {
        for(auto c: m.cells()) xvals[c.id()] = vf[c][d];
        status = TECDAT112( &num_elem, xvals.data(), &VIsDouble );
        assert( status == 0 && "error with TECDAT" );
      } // for
    } // for

    //--------------------------------------------------------------------------
    // WRITE CONNECTIVITY
    //--------------------------------------------------------------------------

    tec_int_t * FaceNodeCounts            = nullptr; // This is NULL for polygonal zones
    tec_int_t * FaceBndryConnectionCounts = nullptr; // not using boundary connections
    tec_int_t * FaceBndryConnectionElems  = nullptr; // not using boundary connections
    tec_int_t * FaceBndryConnectionZones  = nullptr; // not using boundary connections

    // element definitions
    vector<tec_int_t> face_nodes( 2 * num_edges );
    vector<tec_int_t> face_cell_right( num_edges );
    vector<tec_int_t> face_cell_left ( num_edges );

    auto f = 0, i = 0;
    for (auto e : m.edges()) {
      // the nodes of each face
      for (auto v : m.vertices(e))
        face_nodes[i++] = v.id() + 1; // 1-based ids      
      // the cells of each face (1-based ids)
      auto cells = m.cells( e );
      assert( cells.size() > 0 && "cell has no edges");
      // always has left cell
      face_cell_left[f] = cells[0].id() + 1;
      // boundary faces don't have right cell
      face_cell_right[f] = ( cells.size() == 2 ) ? cells[1].id()+1 : 0;
      // incrememnt 
      f++;
    }

    status = TECPOLY112( 
      FaceNodeCounts,
      face_nodes.data(),
      face_cell_left.data(),
      face_cell_right.data(),
      FaceBndryConnectionCounts,
      FaceBndryConnectionElems,
      FaceBndryConnectionZones
    );
    assert( status == 0 && "error with TECNOD" );


    //--------------------------------------------------------------------------
    // CLOSE FILE
    //--------------------------------------------------------------------------

    status = TECEND112();
    assert( status == 0 && "error with TECEND" );

  
    return status;


#else

    std::cerr << "FLECSI not build with tecio support." << std::endl;
    std::exit(1);

    return -1;

#endif

  } // io_tecplot_binary_t::write

    //============================================================================
    //! Implementation of tecplot mesh read for burton specialization.
    //!
    //! \param[in] name Read burton mesh \e m to \e name.
    //! \param[in] m Burton mesh to Read to \e name.
    //!
    //! \return tecplot error code. 0 on success.
    //!
    //============================================================================
  int32_t read( const std::string &name, burton_mesh_t &m)  override
  {
    raise_implemented_error( "No tecplot read functionality has been implemented" );
  };


}; // io_tecplot_ascii_t


////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_tecplot_ascii_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_tecplot_ascii_t.
//!
//! \return Pointer to io_base_t base class of io_tecplot_ascii_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_t> * create_io_tecplot_ascii()
{
  return new burton_io_tecplot_ascii_t;
} // create_io_tecplot_ascii

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_tecplot_binary_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_tecplot_binary_t.
//!
//! \return Pointer to io_base_t base class of io_tecplot_binary_t.
////////////////////////////////////////////////////////////////////////////////
inline flecsi::io_base_t<burton_mesh_t> * create_io_tecplot_binary()
{
  return new burton_io_tecplot_binary_t;
} // create_io_tecplot_binary


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "plt" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_tecplot_dat_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType(
    "plt", create_io_tecplot_binary );

////////////////////////////////////////////////////////////////////////////////
//! Register file extension "dat" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_tecplot_plt_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType(
    "dat", create_io_tecplot_ascii );


} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
