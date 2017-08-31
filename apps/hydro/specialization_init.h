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
#include "../common/specialization_init.h"
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
  // The Closures of other entities
  //----------------------------------------------------------------------------
  
  for ( int i=0; i<num_dims; ++i ) {
    apps::common::color_entity( 
      communicator.get(), 
      closure, 
      remote_info_map, 
      shared_cells_map,
      closure_intersection_map,
      mesh_def.entities(num_dims, i),
      mesh_def.entities(i, num_dims),
      entities[i], 
      entity_color_info[i]
    );
  }


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
      index_spaces::entity_map[0][i], entities[i], coloring_info
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
      ai.from_index_space = index_spaces::entity_map[0][ from_dim ];
      ai.to_index_space = index_spaces::entity_map[0][to_dim ];
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
  apps::common::create_cells( mesh_def, mesh );
 
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
flecsi_register_task(initialize_mesh, loc, single|flecsi::leaf);
  

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver
///////////////////////////////////////////////////////////////////////////////
int specialization_tlt_init(int argc, char ** argv)
{
  
  clog(info) << "In specialization top-level-task init" << std::endl;
  
  // set exceptions 
  apps::common::enable_exceptions();
  
  // get the color
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

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
  if ( !mesh_filename_string.empty() ) {
    if ( rank == 0 )
      std::cout << "Using mesh file \"" << mesh_filename_string << "\"." 
                << std::endl;
  }
  else {
    raise_runtime_error( "No mesh file provided" );
  }
  
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

  return 0;
}


} // namespace
} // namespace
