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
#include "ale/mesh/burton/burton_mesh.h"
#include "ale/utils/errors.h"
#include "ale/utils/string_utils.h"


namespace ale {
namespace mesh {



////////////////////////////////////////////////////////////////////////////////
/// \brief provides base functionality for tecplot writer
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_io_tecplot_base {

public:

  //============================================================================
  //! Typedefs
  //============================================================================

  //! the mesh type
  using mesh_t = burton_mesh_t<N>;

  //! other useful types
  using    size_t = typename mesh_t::size_t;
  using integer_t = typename mesh_t::integer_t;
  using    real_t = typename mesh_t::real_t;
  using   point_t = typename mesh_t::point_t;
  using  vector_t = typename mesh_t::vector_t;
  using  vertex_t = typename mesh_t::vertex_t;
  using    face_t = typename mesh_t::face_t;
  using    cell_t = typename mesh_t::cell_t;

  //! some tecplot typedefs
  using tec_real_t = real_t;
  using tec_int_t = INTEGER4;


  //! class for variable locations
  enum class tec_var_location_t 
  {
    cell = 0,
    node = 1
  };

  //============================================================================
  //! A mapping object utility for tecplot zones
  //============================================================================
  struct tec_zone_map_t {

    //--------------------------------------------------------------------------
    //! \brief constructor
    //--------------------------------------------------------------------------
    tec_zone_map_t( mesh_t & m ) : 
      mesh(m),
      num_zones( m.num_regions() ), 
      elem_zone_map( num_zones )
    {
      std::vector< size_t > local_elem_id( num_zones, 0 );

      // determine a local cell zone ordering
      for ( auto c : m.cells() ) {
        auto reg_id = c->region();
        elem_zone_map[reg_id][c] = local_elem_id[reg_id]++;
      }
    }
      

    //--------------------------------------------------------------------------
    //! \brief create a region map
    //--------------------------------------------------------------------------
    template< typename T >
    void build_face_map( size_t zone_id, const T & elem_this_zone ) {

      // count how many faces there are in this block
      num_faces_this_zone = 0;
      for ( auto c : elem_this_zone ) 
        num_faces_this_zone += mesh.faces(c).size();
      
      // create a face map
      std::vector< face_t* > faces_this_zone; 
      faces_this_zone.reserve( num_faces_this_zone );
      
      for ( auto c : elem_this_zone ) 
        for ( auto f : mesh.faces( c ) )
          faces_this_zone.emplace_back( f );
      
      // sort, then delete duplicate entries
      std::sort( faces_this_zone.begin(), faces_this_zone.end() );
      auto last = std::unique( faces_this_zone.begin(), faces_this_zone.end() );
      faces_this_zone.erase( last, faces_this_zone.end() ); 
      
      // readjust the number of faces
      num_faces_this_zone = faces_this_zone.size();

      // count total face nodes
      size_t num_nodes_this_zone = 0;
      for ( auto f : faces_this_zone ) 
          num_nodes_this_zone += mesh.vertices(f).size();
      
      // face definitions
      face_nodes.resize( num_nodes_this_zone );
      face_cell_right.resize( num_faces_this_zone );
      face_cell_left .resize( num_faces_this_zone );
      
      // for zone-connectivity
      face_conn_counts.clear();
      face_conn_elems.clear();
      face_conn_zones.clear();

      std::cout << "building for " << zone_id << std::endl; 
      
      size_t f = 0, i = 0;
      num_face_conn = 0;
      for (auto face : faces_this_zone) {
        // the nodes of each face
        for (auto vert : mesh.vertices(face))
          face_nodes[i++] = vert.id() + 1; // 1-based ids      
        // the cells of each face (1-based ids)
        auto cells = mesh.cells( face );
        assert( cells.size() > 0 && "cell has no edges");
        // initialize cells to no connections
        face_cell_left [f] = 0;
        face_cell_right[f] = 0;
        // this region id
        auto this_region = zone_id;
        // always has left cell
        auto left_cell = cells[0];
        auto left_region = left_cell->region();
        // left cell is local
        if ( left_region == zone_id )
          face_cell_left[f] = elem_zone_map[ this_region ].at( left_cell ) + 1;
        // left cell is on another zone
        else {
          face_cell_left[f] = - (++num_face_conn);
          face_conn_counts.emplace_back( 1 );
          face_conn_elems.emplace_back( 
            elem_zone_map[ left_region ].at( left_cell ) + 1 );
          face_conn_zones.emplace_back( left_region + 1 );
        }         
        // boundary faces don't have right cell
        if ( cells.size() > 1 ) {
          auto right_cell = cells[1];
          auto right_region = right_cell->region();
          // right cell is local
          if ( right_region == zone_id ) 
            face_cell_right[f] = elem_zone_map[ this_region ].at( right_cell ) + 1;
          // right cell is on another zone
          else {
            face_cell_right[f] = - (++num_face_conn);
            face_conn_counts.emplace_back( 1 );
            face_conn_elems.emplace_back( 
              elem_zone_map[ right_region ].at( right_cell ) + 1 );
            face_conn_zones.emplace_back( right_region );
          }
        }
        // incrememnt 
        f++;
      }
      
    }

    //--------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------

    //! \brief a referemce for the mesh
    mesh_t & mesh;

    //! \brief number of zones
    tec_int_t num_zones;

    //! \brief for face-zone connectivity
    tec_int_t num_faces_this_zone = 0;
    std::vector<tec_int_t> face_nodes;
    std::vector<tec_int_t> face_cell_right;
    std::vector<tec_int_t> face_cell_left;
    
    //! \brief  for elememnt-zone connectivity
    tec_int_t num_face_conn = 0;
    std::vector<tec_int_t> face_conn_counts;
    std::vector<tec_int_t> face_conn_elems;
    std::vector<tec_int_t> face_conn_zones;
    
    //! \brief  storage for the zone-element mapping
    std::vector< 
      std::map< cell_t*, size_t > 
    > elem_zone_map;
        
  };



  //============================================================================
  //! \brief extract the variable information from the 
  //============================================================================
  // auto 
};


////////////////////////////////////////////////////////////////////////////////
/// \class io_tecplot_ascii_t io_tecplot.h
/// \brief io_tecplot_ascii_t provides a derived type of io_base.h and registrations
///   of the tecplot_ascii file extensions.
///
/// \tparam mesh_t Mesh to template io_base_t on.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_io_tecplot_ascii_t : 
    public flecsi::io_base_t< burton_mesh_t<N> >, burton_io_tecplot_base<N>
{

  //============================================================================
  //! typedefs
  //============================================================================
  using base_t = burton_io_tecplot_base<N>;

  using typename base_t::tec_var_location_t;
  using typename base_t::tec_zone_map_t;

  using typename base_t::tec_real_t;
  using typename base_t::tec_int_t;

  using typename base_t::real_t;
  using typename base_t::integer_t;
  using typename base_t::vector_t;

  using typename base_t::mesh_t;
  
  //============================================================================
  //! Default constructor
  //============================================================================
  burton_io_tecplot_ascii_t() = default;


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
  int32_t write( const std::string &name, mesh_t &m) override
  {

    // some aliazes
    using std::endl;
    using std::pair;
    using std::make_pair;
    using std::string;
    using std::vector;


    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return std::forward<decltype(str)>(str);
    };

    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

    std::cout << "Writing mesh to: " << name << std::endl;

    // open the file for writing
    std::ofstream ofs( name.c_str() );
    assert( ofs.good() && "error opening file" );

    //--------------------------------------------------------------------------
    // collect field data
    //--------------------------------------------------------------------------

    // get the general statistics
    auto num_dims  = m.num_dimensions();

    // variable extension for vectors
    std::string var_ext[3];
    var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

    // collect the variables
    vector< pair<string,tec_var_location_t> > variables;

    // write the coordinate names
    if ( num_dims > 0 )
      variables.emplace_back( make_pair( "x", tec_var_location_t::node ) );
    if ( num_dims > 1 )
      variables.emplace_back( make_pair( "y", tec_var_location_t::node ) );
    if ( num_dims > 2 )
      variables.emplace_back( make_pair( "z", tec_var_location_t::node ) );

    //----------------------------------------------------------------------------
    // nodal field data
    //----------------------------------------------------------------------------

    int num_nf = 0; // number of nodal fields
    // real scalars persistent at vertices
    auto rspav = access_type_if(m, real_t, is_persistent_at(m,vertices));
    num_nf += rspav.size();
    // int scalars persistent at vertices
    auto ispav = access_type_if(m, integer_t, is_persistent_at(m,vertices));
    num_nf += ispav.size();
    // real vectors persistent at vertices
    auto rvpav = access_type_if(m, vector_t, is_persistent_at(m,vertices));
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
    auto rspac = access_type_if(m, real_t, is_persistent_at(m,cells));
    num_ef += rspac.size();
    // int scalars persistent at cells
    auto ispac = access_type_if(m, integer_t, is_persistent_at(m,cells));
    num_ef += ispac.size();
    // real vectors persistent at cells
    auto rvpac = access_type_if(m, vector_t, is_persistent_at(m,cells));
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

    //--------------------------------------------------------------------------
    // Element-Zone connectivity
    //--------------------------------------------------------------------------

    auto num_zones = m.num_regions();
   
    // get the regions
    auto region_cells = m.regions();
    
    // create a region map
    tec_zone_map_t mapping( m );

    //--------------------------------------------------------------------------
    // Loop over Regions
    //--------------------------------------------------------------------------
  
    for ( size_t izn=0; izn<num_zones; izn++ ) {

      // get the elements in this block
      const auto & elem_this_zone = region_cells[izn];
      auto num_elem_this_zone = elem_this_zone.size();

      //------------------------------------------------------------------------
      // face/edge connectivity

      mapping.build_face_map( izn, elem_this_zone );
      
      //------------------------------------------------------------------------
      // Zone Header

      ofs << "ZONE" << endl;
      ofs << " T = \"region " << izn << "\"" << endl;
      ofs << " ZONETYPE=FEPOLYGON" << endl;
      ofs << " NODES=" << m.num_vertices() 
          << " FACES=" << mapping.num_faces_this_zone
          << " ELEMENTS=" << num_elem_this_zone << endl;
      ofs << " NumConnectedBoundaryFaces=" << mapping.num_face_conn
          << " TotalNumBoundaryConnections=" << mapping.num_face_conn << endl;
      ofs << " DATAPACKING=BLOCK" << endl;
      
      auto nf_start = 1;
      auto nf_end = num_dims + num_nf;
      auto ef_start = nf_end + 1;
      auto ef_end = nf_end + num_ef;
      if ( num_ef > 0 )
        ofs << " VARLOCATION=([" << ef_start << "-" << ef_end << "]=CELLCENTERED)" 
            << endl;
      
      //------------------------------------------------------------------------
      // Write Nodal Data for Zone 1 only
      if ( izn == 0 ) {
        
        //----------------------------------------------------------------------
        // WRITE coordinates
      
        // get the coordinates from the mesh.
        for(int d=0; d < num_dims; ++d) {
          for (auto v : m.vertices()) 
            ofs << v->coordinates()[d] << " ";
          ofs << endl;
        } // for
        
        //------------------------------------------------------------------------
        // nodal field data
        
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

      }
      //------------------------------------------------------------------------
      // all other zones share nodal data
      else {

        ofs << " VARSHARELIST=([" << nf_start << "-" << nf_end << "]=1)" << endl;         
          
      } // first zone


      //------------------------------------------------------------------------
      // cell field data

      // element field buffer
      for(auto sf: rspac) {
        for(auto c: elem_this_zone) ofs << sf[c] << " ";
        ofs << endl;
      } // for
      for(auto sf: ispac) {
        for(auto c: elem_this_zone) ofs << sf[c] << " ";
        ofs << endl;
      } // for
      for(auto vf: rvpac) {
        for(int d=0; d < num_dims; ++d) {
          for(auto c: elem_this_zone) ofs << vf[c][d] << " ";
          ofs << endl;
        } // for
      } // for

      //------------------------------------------------------------------------
      // WRITE CONNECTIVITY

      ofs << "#face nodes" << endl;
      for (auto i : mapping.face_nodes) ofs << i << " ";
      ofs << endl;
      
      ofs << "#left elements" << endl;
      for (auto i : mapping.face_cell_left) ofs << i << " ";
      ofs << endl;
      
      ofs << "#right elements" << endl;
      for (auto i : mapping.face_cell_right) ofs << i << " ";
      ofs << endl;

      ofs << "#boundary connection counts" << endl;
      for (auto i : mapping.face_conn_counts) ofs << i << " ";
      ofs << endl;

      ofs << "#boundary connection elements" << endl;
      for (auto i : mapping.face_conn_elems) ofs << i << " ";
      ofs << endl;

      ofs << "#boundary connection zones" << endl;
      for (auto i : mapping.face_conn_zones) ofs << i << " ";
      ofs << endl;
 
    } // block

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
  int32_t read( const std::string &name, mesh_t &m)  override
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
template<std::size_t N>
struct burton_io_tecplot_binary_t : 
    public flecsi::io_base_t< burton_mesh_t<N> >, burton_io_tecplot_base<N>
{

  //============================================================================
  //! typedefs
  //============================================================================
  using base_t = burton_io_tecplot_base<N>;

  using typename base_t::tec_var_location_t;
  using typename base_t::tec_zone_map_t;

  using typename base_t::tec_real_t;
  using typename base_t::tec_int_t;

  using typename base_t::real_t;
  using typename base_t::integer_t;
  using typename base_t::vector_t;

  using typename base_t::mesh_t;
  
  //============================================================================
  //! Default constructor
  //============================================================================
  burton_io_tecplot_binary_t() = default;

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
  int32_t write( const std::string &name, mesh_t &m)  override
  {

#ifdef HAVE_TECIO

    std::cout << "Writing mesh to: " << name << std::endl;

    // some aliazes
    using std::endl;
    using std::pair;
    using std::make_pair;
    using std::string;
    using std::vector;

    // a lambda function for validating strings
    auto validate_string = []( auto && str ) {
      return utils::replace_all( std::forward<decltype(str)>(str), " ", "_" );;
    };


    //--------------------------------------------------------------------------
    // setup
    //--------------------------------------------------------------------------

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
    tec_int_t num_elem  = m.num_cells();
    tec_int_t num_zones = m.num_regions();

    // set the time
    double soln_time = m.time();



    //--------------------------------------------------------------------------
    // collect field data
    //--------------------------------------------------------------------------

    // variable extension for vectors
    std::string var_ext[3];
    var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

    // collect the variables
    vector< pair<string,tec_var_location_t> > variables;

    // write the coordinate names
    if ( num_dims > 0 )
      variables.emplace_back( make_pair( "x", tec_var_location_t::node ) );
    if ( num_dims > 1 )
      variables.emplace_back( make_pair( "y", tec_var_location_t::node ) );
    if ( num_dims > 2 )
      variables.emplace_back( make_pair( "z", tec_var_location_t::node ) );


    //----------------------------------------------------------------------------
    // nodal field data

    // real scalars persistent at vertices
    auto rspav = access_type_if(m, real_t, is_persistent_at(m,vertices));
    // int scalars persistent at vertices
    auto ispav = access_type_if(m, integer_t, is_persistent_at(m,vertices));
    // real vectors persistent at vertices
    auto rvpav = access_type_if(m, vector_t, is_persistent_at(m,vertices));

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
    auto rspac = access_type_if(m, real_t, is_persistent_at(m,cells));
    // int scalars persistent at cells
    auto ispac = access_type_if(m, integer_t, is_persistent_at(m,cells));
    // real vectors persistent at cells
    auto rvpac = access_type_if(m, vector_t, is_persistent_at(m,cells));


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




    //--------------------------------------------------------------------------
    // Var locations
    //--------------------------------------------------------------------------

    // create the variable location array.  also, we write all nodal
    // data to the first zone, and share the rest
    vector<tec_int_t> var_locations, var_sharing;
    for ( const auto & var : variables ) {
      auto var_loc = var.second;
      // cast var location to an integer
      var_locations.push_back( static_cast<tec_int_t>(var_loc) );
      // 0 means not shared, otherwise, share nodal data with first zone
      auto share_var =  (var_loc == tec_var_location_t::node) ? 1 : 0;
      var_sharing.push_back( share_var );
    }

    //--------------------------------------------------------------------------
    // Element-Zone connectivity
    //--------------------------------------------------------------------------

    // get the regions
    auto region_cells = m.regions();

    // create a region map
    tec_zone_map_t mapping( m );

    //--------------------------------------------------------------------------
    // Loop over Regions
    //--------------------------------------------------------------------------
  
    for ( size_t izn=0; izn<num_zones; izn++ ) {

      // get the elements in this block
      const auto & elem_this_zone = region_cells[izn];
      tec_int_t num_elem_this_zone = elem_this_zone.size();

      //------------------------------------------------------------------------
      // face/edge connectivity
      
      mapping.build_face_map( izn, elem_this_zone );

      //------------------------------------------------------------------------
      // Create ZONE header

      std::string ZoneTitle = "region " + std::to_string(izn);

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
      tec_int_t TotalNumFaceNodes        = 2*mapping.num_faces_this_zone; // total nodes for all faces
      tec_int_t TotalNumBndryFaces       = mapping.num_face_conn; // number boundary connections
      tec_int_t TotalNumBndryConnections = mapping.num_face_conn; // number boundary connections
      tec_int_t ShareConnectivityFromZone = 0;  // pass 0 to indicate no sharign
      tec_int_t * PassiveVarList = nullptr;
      tec_int_t * ShareVarFromZone = nullptr;

      // share nodal data from first zone
      if ( izn > 0 ) ShareVarFromZone = var_sharing.data();

      status = TECZNE112( const_cast<char*>( ZoneTitle.c_str() ),
                          &ZoneType,
                          &num_nodes,
                          &num_elem_this_zone,
                          &mapping.num_faces_this_zone,
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



      // a temporary for vector values
      vector<tec_real_t> xvals( std::max( num_nodes, num_elem ) );
      vector<tec_real_t> yvals( std::max( num_nodes, num_elem ) );
      

      //------------------------------------------------------------------------
      // Write Nodal Data for Zone 1 only
      if ( izn == 0 ) {
        
        //----------------------------------------------------------------------
        // write coordinates
        
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

        //----------------------------------------------------------------------
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

      } // first zone


      //------------------------------------------------------------------------
      // cell field data

      // element field buffer
      for(auto sf: rspac) {
        size_t cid = 0;
        for(auto c: elem_this_zone) xvals[cid++] = sf[c];
        status = TECDAT112( &num_elem_this_zone, xvals.data(), &VIsDouble );
        assert( status == 0 && "error with TECDAT" );
      } // for
      for(auto sf: ispac) {
        // cast int fields to real_t
        size_t cid = 0;
        for(auto c: elem_this_zone) xvals[cid++] = (tec_real_t)sf[c];
        status = TECDAT112( &num_elem_this_zone, xvals.data(), &VIsDouble );
        assert( status == 0 && "error with TECDAT" );
      } // for
      for(auto vf: rvpac) {
        for(int d=0; d < num_dims; ++d) {
          size_t cid = 0;
          for(auto c: elem_this_zone) xvals[cid++] = vf[c][d];
          status = TECDAT112( &num_elem_this_zone, xvals.data(), &VIsDouble );
          assert( status == 0 && "error with TECDAT" );
        } // for
      } // for

      //--------------------------------------------------------------------------
      // WRITE CONNECTIVITY

      tec_int_t * FaceNodeCounts = nullptr; // This is NULL for polygonal zones

      status = TECPOLY112( 
        FaceNodeCounts,
        mapping.face_nodes.data(),
        mapping.face_cell_left.data(),
        mapping.face_cell_right.data(),
        mapping.face_conn_counts.data(),
        mapping.face_conn_elems.data(),
        mapping.face_conn_zones.data()
      );
      assert( status == 0 && "error with TECNOD" );


    } // blocks

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
  int32_t read( const std::string &name, mesh_t &m)  override
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
template< std::size_t N >
inline flecsi::io_base_t< burton_mesh_t<N> > * create_io_tecplot_ascii()
{
  return new burton_io_tecplot_ascii_t<N>;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Create an io_tecplot_binary_t and return a pointer to the base class.
//!
//! \tparam mesh_t Mesh type for io_tecplot_binary_t.
//!
//! \return Pointer to io_base_t base class of io_tecplot_binary_t.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
inline flecsi::io_base_t< burton_mesh_t<N> > * create_io_tecplot_binary()
{
  return new burton_io_tecplot_binary_t<N>;
}


////////////////////////////////////////////////////////////////////////////////
//! Register file extension "dat" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_2d_tecplot_dat_registered =
  flecsi::io_factory_t< burton_mesh_t<2> >::instance().registerType(
    "dat", create_io_tecplot_ascii<2> );

static bool burton_3d_tecplot_dat_registered =
  flecsi::io_factory_t< burton_mesh_t<3> >::instance().registerType(
    "dat", create_io_tecplot_ascii<3> );

////////////////////////////////////////////////////////////////////////////////
//! Register file extension "plt" with factory.
////////////////////////////////////////////////////////////////////////////////
static bool burton_2d_tecplot_plt_registered =
  flecsi::io_factory_t< burton_mesh_t<2> >::instance().registerType(
    "plt", create_io_tecplot_binary<2> );

static bool burton_3d_tecplot_plt_registered =
  flecsi::io_factory_t< burton_mesh_t<3> >::instance().registerType(
    "plt", create_io_tecplot_binary<3> );


} // namespace mesh
} // namespace ale

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
