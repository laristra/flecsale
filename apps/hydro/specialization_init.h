/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This is the main initialization driver
///////////////////////////////////////////////////////////////////////////////
#pragma once

// user incldues
#include "../common/exceptions.h"
#include "arguments.h"
#include "cinch/logging/cinchlog.h"
#include "flecsi/coloring/dcrs_utils.h"
#include "flecsi/coloring/mpi_communicator.h"
#include "flecsi/coloring/parmetis_colorer.h"
#include "flecsi/execution/execution.h"
#include "flecsi/io/exodus_definition.h"

// system includes
#include <iostream>

namespace apps {
namespace hydro {
  
////////////////////////////////////////////////////////////////////////////////
//! Build the list of exclusve, shared and ghost entities
////////////////////////////////////////////////////////////////////////////////
template < 
  int entity_dim,
  typename MESH_DEF, 
  typename COMM,
  typename CLOSURE_SET,
  typename ENTITY_MAP,
  typename INTERSECTION_MAP,
  typename INDEX_COLOR,
  typename COLOR_INFO,
  int DIM = MESH_DEF::dimension()
>
auto color_entity(
  const MESH_DEF & mesh_def, 
  COMM * communicator,
  const CLOSURE_SET & closure, 
  const ENTITY_MAP & remote_info_map,
  const ENTITY_MAP & shared_cells_map,
  const INTERSECTION_MAP & closure_intersection_map,
  INDEX_COLOR & entities, 
  COLOR_INFO & entity_color_info
) {

  // some compile time constants
  constexpr auto cell_dim = DIM;
  
  // some type aliases
  using entity_info_t = flecsi::coloring::entity_info_t;

  // info about the mpi communicator
  auto comm_size = communicator->size();
  auto rank = communicator->rank();

  // Form the entity closure
  auto entity_closure = 
    flecsi::topology::entity_closure<cell_dim, entity_dim>(mesh_def, closure);

  // Assign entity ownership
  std::vector<std::set<size_t>> entity_requests(comm_size);
  std::set<entity_info_t> entity_info;

  {
    size_t offset(0);
    for(auto i: entity_closure) {

      // Get the set of cells that reference this entity.
      auto referencers = 
        flecsi::topology::entity_referencers<cell_dim, entity_dim>(mesh_def, i);

      size_t min_rank(std::numeric_limits<size_t>::max());
      std::set<size_t> shared_entities;

      // Iterate the direct referencers to assign entity ownership.
      for(auto c: referencers) {

        // Check the remote info map to see if this cell is
        // off-color. If it is, compare it's rank for
        // the ownership logic below.
        if(remote_info_map.find(c) != remote_info_map.end()) {
          min_rank = std::min(min_rank, remote_info_map.at(c).rank);
          shared_entities.insert(remote_info_map.at(c).rank);
        }
        else {
          // If the referencing cell isn't in the remote info map
          // it is a local cell.

          // Add our rank to compare for ownership.
          min_rank = std::min(min_rank, size_t(rank));

          // If the local cell is shared, we need to add all of
          // the ranks that reference it.
          if(shared_cells_map.find(c) != shared_cells_map.end()) 
            shared_entities.insert( 
              shared_cells_map.at(c).shared.begin(),
              shared_cells_map.at(c).shared.end()
            );
        } // if

        // Iterate through the closure intersection map to see if the
        // indirect reference is part of another rank's closure, i.e.,
        // that it is an indirect dependency.
        for(auto ci: closure_intersection_map) 
          if(ci.second.find(c) != ci.second.end()) 
            shared_entities.insert(ci.first);
      } // for

      if(min_rank == rank) {
        // This is a entity that belongs to our rank.
        auto entry = entity_info_t(i, rank, offset++, shared_entities);
        entity_info.insert( entry );
      }
      else {
        // Add remote entity to the request for offset information.
        entity_requests[min_rank].insert(i);
      } // if
    } // for
  } // scope

  auto entity_offset_info =
    communicator->get_entity_info(entity_info, entity_requests);

  // Vertices index coloring.
  for(auto i: entity_info) {
    // if it belongs to other colors, its a shared entity
    if(i.shared.size()) {
      entities.shared.insert(i);
      // Collect all colors with whom we require communication
      // to send shared information.
      entity_color_info.shared_users = 
        flecsi::utils::set_union(entity_color_info.shared_users, i.shared);
    }
    // otherwise, its exclusive
    else 
      entities.exclusive.insert(i);
  } // for

  {
    size_t r(0);
    for(auto i: entity_requests) {

      auto offset(entity_offset_info[r].begin());
      for(auto s: i) {
        entities.ghost.insert(entity_info_t(s, r, *offset));
        // Collect all colors with whom we require communication
        // to receive ghost information.
        entity_color_info.ghost_owners.insert(r);
        // increment counter
        ++offset;
      } // for

      ++r;
    } // for
  } // scope
  
  entity_color_info.exclusive = entities.exclusive.size();
  entity_color_info.shared = entities.shared.size();
  entity_color_info.ghost = entities.ghost.size();
}

////////////////////////////////////////////////////////////////////////////////
/// \brief the main cell coloring driver
////////////////////////////////////////////////////////////////////////////////
void partition_mesh( char_array_t filename ) 
{
  // set some compile time constants 
  constexpr auto num_dims = mesh_t::num_dimensions;
  constexpr auto cell_dim = num_dims;
  constexpr auto thru_dim = 0;

  // make some type aliases
  using real_t = mesh_t::real_t;
  using exodus_definition_t = flecsi::io::exodus_definition__<num_dims, real_t>;
  using entity_info_t = flecsi::coloring::entity_info_t;
  
  // load the mesh
  auto filename_string = filename.str();
  exodus_definition_t mesh_def( filename_string );

  // Create a communicator instance to get neighbor information.
  auto communicator = std::make_unique<flecsi::coloring::mpi_communicator_t>();
  auto comm_size = communicator->size();
  auto rank = communicator->rank();


  // create a vector of colorings and color info for each dimensional entity
  std::vector< flecsi::coloring::index_coloring_t > 
    entities( num_dims+1 );
  std::vector< flecsi::coloring::coloring_info_t >
    entity_color_info( num_dims+1 );

  //----------------------------------------------------------------------------
  // Cell Coloring
  //----------------------------------------------------------------------------
    
  // Cells index coloring.
  auto & cells = entities[ num_dims ];
  auto & cell_color_info = entity_color_info[ num_dims ];
 
  // Create the dCRS representation for the distributed colorer.
  // This essentialy makes the graph of the dual mesh.
  auto dcrs = flecsi::coloring::make_dcrs(mesh_def);

  // Create a colorer instance to generate the primary coloring.
  auto colorer = std::make_unique<flecsi::coloring::parmetis_colorer_t>();

  // Create the primary coloring.
  cells.primary = colorer->color(dcrs);

  //----------------------------------------------------------------------------
  // Cell Closure.  However many layers of ghost cells are needed are found 
  // here.
  //----------------------------------------------------------------------------
    
  // Compute the dependency closure of the primary cell coloring
  // through vertex intersections (specified by last argument "1").
  // To specify edge or face intersections, use 2 (edges) or 3 (faces).
  auto closure = 
    flecsi::topology::entity_neighbors<cell_dim, cell_dim, thru_dim>(
      mesh_def, cells.primary
    );

  // Subtracting out the initial set leaves just the nearest
  // neighbors. This is similar to the image of the adjacency
  // graph of the initial indices.
  auto nearest_neighbors = 
    flecsi::utils::set_difference(closure, cells.primary);

  //----------------------------------------------------------------------------
  // Find one more layer of ghost cells now, these are needed to get some corner
  // cases right.
  //----------------------------------------------------------------------------

  // We can iteratively add halos of nearest neighbors, e.g.,
  // here we add the next nearest neighbors. For most mesh types
  // we actually need information about the ownership of these indices
  // so that we can deterministically assign rank ownership to vertices.
  auto nearest_neighbor_closure =
    flecsi::topology::entity_neighbors<cell_dim, cell_dim, thru_dim>(
      mesh_def, nearest_neighbors
    );

  // Subtracting out the closure leaves just the
  // next nearest neighbors.
  auto next_nearest_neighbors =
    flecsi::utils::set_difference(nearest_neighbor_closure, closure);

  // The union of the nearest and next-nearest neighbors gives us all
  // of the cells that might reference a vertex that we need.
  auto all_neighbors = 
    flecsi::utils::set_union(nearest_neighbors, next_nearest_neighbors);

  //----------------------------------------------------------------------------
  // Find exclusive, shared, and ghost cells..
  //----------------------------------------------------------------------------

  // Get the rank and offset information for our nearest neighbor
  // dependencies. This also gives information about the ranks
  // that access our shared cells.
  auto cell_nn_info =
    communicator->get_primary_info(cells.primary, nearest_neighbors);

  // Create a map version of the local info for lookups below.
  std::unordered_map<size_t, size_t> primary_indices_map;
  {
    size_t offset(0);
    for(auto i: cells.primary) 
      primary_indices_map[offset++] = i;
  } // scope

  // Populate exclusive and shared cell information.
  {
    size_t offset(0);
    for(auto i: std::get<0>(cell_nn_info)) {
      // get the entity info
      auto entry =
        entity_info_t(primary_indices_map.at(offset), rank, offset, i);
      // if it belongs to other colors, its a shared cell
      if(i.size()) {
        cells.shared.insert(entry);
        // Collect all colors with whom we require communication
        // to send shared information.
        cell_color_info.shared_users = 
          flecsi::utils::set_union(cell_color_info.shared_users, i);
      }
      // otherwise, its an exclusive cell
      else 
        cells.exclusive.insert(entry);
      // increment counter
      offset++;
    } // for
  } // scope
  
  // Populate ghost cell information.
  for(auto i: std::get<1>(cell_nn_info)) {
    cells.ghost.insert(i);
    // Collect all colors with whom we require communication
    // to receive ghost information.
    cell_color_info.ghost_owners.insert(i.rank);
  }
  
  // store the sizes of each set
  cell_color_info.exclusive = cells.exclusive.size();
  cell_color_info.shared = cells.shared.size();
  cell_color_info.ghost = cells.ghost.size();

  //----------------------------------------------------------------------------
  // Create some maps for easy lookups when determining the other dependent
  // closures.
  //----------------------------------------------------------------------------

  // Create a map version for lookups below.
  std::unordered_map<size_t, entity_info_t> shared_cells_map;
  for(auto i: cells.shared)
    shared_cells_map[i.id] = i;

  // Get the rank and offset information for all relevant neighbor
  // dependencies. This information will be necessary for determining
  // shared vertices.
  auto cell_all_info = 
    communicator->get_primary_info(cells.primary, all_neighbors);

  // Create a map version of the remote info for lookups below.
  std::unordered_map<size_t, entity_info_t> remote_info_map;
  for(auto i: std::get<1>(cell_all_info)) 
    remote_info_map[i.id] = i;

  // Get the intersection of our nearest neighbors with the nearest
  // neighbors of other ranks. This map of sets will only be populated
  // with intersections that are non-empty
  auto closure_intersection_map =
    communicator->get_intersection_info(nearest_neighbors);

  //----------------------------------------------------------------------------
  // Vertex Closure
  //----------------------------------------------------------------------------
  
  color_entity<0>( 
    mesh_def, 
    communicator.get(), 
    closure, 
    remote_info_map, 
    shared_cells_map,
    closure_intersection_map,
    entities[0], 
    entity_color_info[0]
  );

  //----------------------------------------------------------------------------
  // Edge Closure
  //----------------------------------------------------------------------------
  
  color_entity<1>( 
    mesh_def, 
    communicator.get(), 
    closure, 
    remote_info_map, 
    shared_cells_map,
    closure_intersection_map,
    entities[1], 
    entity_color_info[1]
  );

  //----------------------------------------------------------------------------
  // Face Closure
  //----------------------------------------------------------------------------
 
  if ( num_dims > 2 )
    color_entity<2>( 
      mesh_def, 
      communicator.get(), 
      closure, 
      remote_info_map, 
      shared_cells_map,
      closure_intersection_map,
      entities[2], 
      entity_color_info[2]
    );


  //----------------------------------------------------------------------------
  // Add the results to the context
  //----------------------------------------------------------------------------

  // Alias the index spaces type
  using index_spaces = mesh_t::index_spaces_t;

  // Get the context instance.
  auto & context = flecsi::execution::context_t::instance();

  // Gather the coloring info from all colors
  for ( int i=0; i<num_dims+1; ++i ) {
    auto coloring_info = 
      communicator->gather_coloring_info(entity_color_info[i]);
    context.add_coloring( 
      index_spaces::entity_map[i], entities[i], coloring_info
    );
  }
    
  //----------------------------------------------------------------------------
  // add adjacency information
  //----------------------------------------------------------------------------
  
  std::vector< std::vector<size_t> > entity_ids(num_dims+1);

  // create a master list of all entities
  for ( int i=0; i<num_dims+1; ++i ) {

    const auto & these_entities = entities[i];
    auto & these_ids = entity_ids[i];

    auto num_entities = 
      these_entities.exclusive.size() + 
      these_entities.shared.size() +
      these_entities.ghost.size();

    these_ids.reserve( num_entities );

    for ( auto e : these_entities.exclusive ) these_ids.push_back(e.id);
    for ( auto e : these_entities.shared    ) these_ids.push_back(e.id);
    for ( auto e : these_entities.ghost     ) these_ids.push_back(e.id);

    // sort the entities ( this is needed for the search down below )
    // hoepefully it doesnt cause any problems
    std::sort( these_ids.begin(), these_ids.end() );
    auto last = std::unique( these_ids.begin(), these_ids.end() );
    if ( last != these_ids.end() )
      clog_error( "Duplicate ids in master lists" );

  }

  // loop over each dimension and determine the adjacency sizes
  for ( int from_dim = 0; from_dim<=num_dims; ++from_dim ) {
   
    // the master list of all entity ids
    const auto & from_ids = entity_ids[from_dim];

    for ( int to_dim = 0; to_dim<=num_dims; ++to_dim ) {

      // skip the case where both dimensions are the same
      if ( from_dim == to_dim ) continue;

      // the master list of all entity ids
      const auto & to_ids = entity_ids[to_dim];
      const auto to_ids_begin = to_ids.begin();
      const auto to_ids_end = to_ids.end();

      // populate the adjacency information
      flecsi::coloring::adjacency_info_t ai;
      ai.index_space = 
        index_spaces::connectivity_map[ from_dim ][ to_dim ];
      ai.from_index_space = index_spaces::entity_map[ from_dim ];
      ai.to_index_space = index_spaces::entity_map[to_dim ];
      ai.color_sizes.resize(comm_size);
  
      // loop over all cells and count the number of adjacencies
      size_t cnt = 0;
      for ( auto c : from_ids ) {
        // get the attached sub entitites
        const auto & ids = mesh_def.entities(from_dim, to_dim, c);
        // we need to make sure they are in this colors master
        // list though
        for ( auto v : ids ) {
          auto it = std::lower_bound( to_ids_begin, to_ids_end, v );
          if ( it != to_ids_end && *it == v ) 
            cnt++;
        }
      }
  
      // gather the results
      ai.color_sizes = communicator->gather_sizes( cnt );
    
      // add the result to the context
      context.add_adjacency(ai);

    }
  }


  //----------------------------------------------------------------------------
  // add intermediate mappings
  //
  // These are needed so that entities that are implicitly created by
  // mesh.init<>() maintain consistent orderings that match the pre-computed
  // colorings.
  //----------------------------------------------------------------------------
  
  for ( int i=1; i<num_dims; ++i ) {

    std::unordered_map<size_t, std::vector<size_t>> entity_to_vertex_map;

    for ( auto e : entity_ids[i] ) { 
      auto vs = mesh_def.entities( /* dimension */ i, /* domain */ 0, e );
      entity_to_vertex_map.emplace( e, vs ); 
    }

    context.add_intermediate_map(
      /* dimension */ i, /* domain */ 0, entity_to_vertex_map
    );

  }

  //----------------------------------------------------------------------------
  // output the result
  //----------------------------------------------------------------------------

  // figure out this ranks file name
  auto basename = utils::basename( filename_string );
  auto output_prefix = utils::remove_extension( basename );
  auto output_filename = output_prefix + "-partition_rank" +
    apps::common::zero_padded(rank) + ".exo";

  // a lamda function to convert sets of entitiy_info_t's to vectors
  // of ids
  auto to_vec = [](const auto & list_in)
  {
    std::vector<size_t> list_out;
    list_out.reserve( list_in.size() );
    for ( auto & e : list_in )
      list_out.push_back( e.id );
    return list_out;
  };

  // get a reference to the vertices
  const auto & vertices = entities[0];

  // open the exodus file
  if ( rank == 0 )
    std::cout << "Writing mesh to: " << output_filename << std::endl;

  using std::make_pair;
  mesh_def.write(
    output_filename,
    { 
      make_pair( "exclusive cells", to_vec(cells.exclusive) ),
      make_pair( "shared cells", to_vec(cells.shared) ),
      make_pair( "ghost cells", to_vec(cells.ghost) )
    },
    { 
      make_pair( "exclusive vertices", to_vec(vertices.exclusive) ),
      make_pair( "shared vertices", to_vec(vertices.shared) ),
      make_pair( "ghost vertices", to_vec(vertices.ghost) )
    }
  );

  clog(info) << "Finished mesh partitioning." << std::endl;


} // somerhing

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 2D version
////////////////////////////////////////////////////////////////////////////////
template<
  int D,
  typename MESH_DEF,
  typename MESH,
  bool Enabled = (D == 2),
  typename = std::enable_if_t< Enabled >
>
void create_cells( MESH_DEF && mesh_def, MESH && mesh )
{
  
  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------
  
  // some constant expressions
  constexpr auto num_dims = mesh_t::num_dimensions;
  
  // Alias the index spaces type
  using index_spaces = mesh_t::index_spaces_t;

  // alias some other types
  using point_t = typename mesh_t::point_t;
  using vertex_t = typename mesh_t::vertex_t;
  using cell_t = typename mesh_t::cell_t;
  
  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // get the entity maps
  // - lid = local id - the id of the entity local to this processor
  // - mid = mesh id - the original id of the entity ( usually from file )
  // - gid = global id - the new global id of the entity ( set by flecsi )
  const auto & vertex_lid_to_mid = context.index_map( index_spaces::vertices );
  const auto & cell_lid_to_mid = context.index_map( index_spaces::cells );
  
  const auto & vertex_mid_to_lid = context.reverse_index_map( 
    index_spaces::vertices
  );

  //----------------------------------------------------------------------------
  // create the vertices
  //----------------------------------------------------------------------------

  std::vector< vertex_t * > vertices;
  vertices.reserve( vertex_lid_to_mid.size() );

  for(auto & vm: vertex_lid_to_mid) { 
    // get the point
    const auto & p = mesh_def.template vertex<point_t>( vm.second );
    // now create it
    auto v = mesh.create_vertex( p );
    vertices.emplace_back(v);
  } // for vertices

  //----------------------------------------------------------------------------
  // create the cells
  //----------------------------------------------------------------------------

  // create the cells
  for(auto & cm: cell_lid_to_mid) { 
    // get the list of vertices
    auto vs = 
      mesh_def.entities( cell_t::dimension, vertex_t::dimension, cm.second );
    // create a list of vertex pointers
    std::vector< vertex_t * > elem_vs( vs.size() );
    // transform the list of vertices to mesh ids
    std::transform(
      vs.begin(),
      vs.end(),
      elem_vs.begin(),
      [&](auto v) { return vertices[ vertex_mid_to_lid.at(v) ]; }
    );
    // create the cell
    auto c = mesh.create_cell( elem_vs ); 
  }
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 3D version
////////////////////////////////////////////////////////////////////////////////
template<
  int D,
  typename MESH_DEF,
  typename MESH,
  bool Enabled = (D == 3),
  typename std::enable_if_t< Enabled >* = nullptr
>
void create_cells( MESH_DEF && mesh_def, MESH && mesh )
{
  
  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------
  
  // some constant expressions
  constexpr auto num_dims = mesh_t::num_dimensions;
  
  // Alias the index spaces type
  using index_spaces = mesh_t::index_spaces_t;

  // alias some other types
  using point_t = typename mesh_t::point_t;
  using vertex_t = typename mesh_t::vertex_t;
  using face_t = typename mesh_t::face_t;
  using cell_t = typename mesh_t::cell_t;
  
  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // get the entity maps
  // - lid = local id - the id of the entity local to this processor
  // - mid = mesh id - the original id of the entity ( usually from file )
  // - gid = global id - the new global id of the entity ( set by flecsi )
  const auto & vertex_lid_to_mid = context.index_map( index_spaces::vertices );
  const auto & face_lid_to_mid = context.index_map( index_spaces::faces );
  const auto & cell_lid_to_mid = context.index_map( index_spaces::cells );

  const auto & vertex_mid_to_lid = 
    context.reverse_index_map( index_spaces::vertices );

  const auto & face_mid_to_lid = 
    context.reverse_index_map( index_spaces::faces );

  

  //----------------------------------------------------------------------------
  // create the vertices
  //----------------------------------------------------------------------------

  std::vector< vertex_t * > vertices;
  vertices.reserve( vertex_lid_to_mid.size() );

  for(auto & vm: vertex_lid_to_mid) { 
    // get the point
    const auto & p = mesh_def.template vertex<point_t>( vm.second );
    // now create it
    auto v = mesh.create_vertex( p );
    vertices.emplace_back(v);
  } // for vertices


  //----------------------------------------------------------------------------
  // create the faces
  //----------------------------------------------------------------------------

  // storage for the face pointers ( not used in 2d )
  std::vector< face_t * > faces;
  
  // reserve size
  faces.reserve( face_lid_to_mid.size() );

  // loop over the faces
  for(auto & fm: face_lid_to_mid) { 
    // get the list of vertices
    auto vs = 
      mesh_def.entities( face_t::dimension, vertex_t::dimension, fm.second );
    // create a list of vertex pointers
    std::vector< vertex_t * > elem_vs( vs.size() );
    // transform the list of vertices to mesh ids
    std::transform(
      vs.begin(),
      vs.end(),
      elem_vs.begin(),
      [&](auto v) { return vertices[ vertex_mid_to_lid.at(v) ]; }
    );
    // create the face
    auto f = mesh.create_face( elem_vs ); 
    faces.emplace_back( f );
  }

  //----------------------------------------------------------------------------
  // create the cells
  //----------------------------------------------------------------------------

  // create the cells
  for(auto & cm: cell_lid_to_mid) { 
    // get the list of faces
    auto fs = 
      mesh_def.entities( cell_t::dimension, face_t::dimension, cm.second );
    // create a list of face pointers
    std::vector< face_t * > elem_fs( fs.size() );
    // transform the list  mesh ids to face pointers
    std::transform(
      fs.begin(),
      fs.end(),
      elem_fs.begin(),
      [&](auto f) { return faces[ face_mid_to_lid.at(f) ]; }
    );
    // create the cell
    auto c = mesh.create_cell( elem_fs ); 
  }
}

////////////////////////////////////////////////////////////////////////////////
/// \brief the main mesh initialization driver
////////////////////////////////////////////////////////////////////////////////
void initialize_mesh( 
  client_handle_w__<mesh_t> mesh, 
  char_array_t filename
) {
  
  //----------------------------------------------------------------------------
  // Fill the mesh information
  //----------------------------------------------------------------------------

  // some constant expressions
  constexpr auto num_dims = mesh_t::num_dimensions;
  
  // alias some types
  using real_t = mesh_t::real_t;
  using exodus_definition_t = flecsi::io::exodus_definition__<num_dims, real_t>;

  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // Load the mesh
  auto filename_string = filename.str();
  exodus_definition_t mesh_def( filename_string );

  // fill the mesh
  create_cells<num_dims>( mesh_def, mesh );
 
  // initialize the mesh
  mesh.init(); 
  if ( rank == 0 )
  mesh.is_valid();

  //----------------------------------------------------------------------------
  // Some debug
  //----------------------------------------------------------------------------
  
  // figure out this ranks file name
  auto basename = utils::basename( filename_string );
  auto output_prefix = utils::remove_extension( basename );
  auto output_filename = output_prefix + "-connectivity_rank" +
    apps::common::zero_padded(rank) + ".txt";

  // dump to file
  if ( rank == 0 )
    std::cout << "Dumping connectivity to: " << output_filename << std::endl;
  std::ofstream file( output_filename );
  mesh.dump( file );

  // close file
  file.close();

}

///////////////////////////////////////////////////////////////////////////////
// Task Registration
///////////////////////////////////////////////////////////////////////////////
flecsi_register_mpi_task(partition_mesh);
flecsi_register_task(initialize_mesh, loc, single);
  

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver
///////////////////////////////////////////////////////////////////////////////
int specialization_tlt_init(int argc, char ** argv)
{
  
  clog(info) << "In specialization top-level-task init" << std::endl;
  
  // set exceptions 
  apps::common::enable_exceptions();

  //===========================================================================
  // Parse arguments
  //===========================================================================
  
  auto args = process_arguments( argc, argv );
  
  // process the simple ones
  if ( args.count("h") )
    return 0;
 
  // get the input file
  auto mesh_filename_string = 
    args.count("m") ? args.at("m") : std::string();

  // override any inputs if need be
  if ( !mesh_filename_string.empty() )
    std::cout << "Using mesh file \"" << mesh_filename_string << "\"." 
              << std::endl;
  else
    raise_runtime_error( "No mesh file provided" );
  
  //===========================================================================
  // Partition mesh
  //===========================================================================
  
  clog(info) << "Partitioning mesh" << std::endl;
  
  // need to put the filename into a statically sized character array
  auto mesh_filename = utils::to_trivial_string( mesh_filename_string );

  // execute the mpi task to partition the mesh
  flecsi_execute_mpi_task(partition_mesh, mesh_filename);
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver
///////////////////////////////////////////////////////////////////////////////
int specialization_spmd_init(int argc, char ** argv)
{
  //===========================================================================
  // Parse arguments
  //===========================================================================
  
  auto args = process_arguments( argc, argv );
  // Assume arguments are sanitized
  
  // get the input file
  auto mesh_filename_string = 
    args.count("m") ? args.at("m") : std::string();
  
  // need to put the filename into a statically sized character array
  auto mesh_filename = utils::to_trivial_string( mesh_filename_string );
  
  //===========================================================================
  // Load the mesh
  //===========================================================================

  // get a mesh handle and call the initialization task
  auto mesh_handle = flecsi_get_client_handle(mesh_t, meshes, mesh0);
  auto f1 = 
    flecsi_execute_task(initialize_mesh, single, mesh_handle, mesh_filename);
  f1.wait();

  return 0;
}


} // namespace
} // namespace
