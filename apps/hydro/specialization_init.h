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

// register a clog tag
clog_register_tag(coloring);

namespace apps {
namespace hydro {
  
////////////////////////////////////////////////////////////////////////////////
//! Build the list of exclusve, shared and ghost entities
////////////////////////////////////////////////////////////////////////////////
template < 
  int entity_dim,
  typename MD, 
  typename COMM,
  typename CLOSURE_SET,
  typename ENTITY_MAP,
  typename INTERSECTION_MAP,
  typename INDEX_COLOR,
  typename COLOR_INFO,
  int DIM = MD::dimension()
>
auto color_entity(
  const MD & md, 
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
    flecsi::topology::entity_closure<cell_dim, entity_dim>(md, closure);

  // Assign entity ownership
  std::vector<std::set<size_t>> entity_requests(comm_size);
  std::set<entity_info_t> entity_info;

  {
    size_t offset(0);
    for(auto i: entity_closure) {

      // Get the set of cells that reference this entity.
      auto referencers = 
        flecsi::topology::entity_referencers<cell_dim, entity_dim>(md, i);

      {
        clog_tag_guard(coloring);
        clog_container_one(info, i << " referencers", referencers, clog::space);
      } // guard

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

  {
    clog_tag_guard(coloring);
    clog_container_one(
      info, 
      "exclusive entities("<<entity_dim<<")", 
      entities.exclusive, 
      clog::newline
    );
    clog_container_one(
      info, "shared entities("<<entity_dim<<")", entities.shared, clog::newline
    );
    clog_container_one(
      info, "ghost entities("<<entity_dim<<")" , entities.ghost, clog::newline
    );
  } // guard

  
  entity_color_info.exclusive = entities.exclusive.size();
  entity_color_info.shared = entities.shared.size();
  entity_color_info.ghost = entities.ghost.size();
}

////////////////////////////////////////////////////////////////////////////////
/// \brief the main cell coloring driver
////////////////////////////////////////////////////////////////////////////////
void partition_mesh( char_array_t filename ) 
{
  clog(info) << "Starting mesh partitioner." << std::endl;

  // set some compile time constants 
  constexpr auto num_dims = mesh_t::num_dimensions;
  constexpr auto cell_dim = num_dims;
  constexpr auto thru_dim = 0;

  // make some type aliases
  using real_t = mesh_t::real_t;
  using exodus_definition_t = flecsi::io::exodus_definition__<num_dims, real_t>;
  using entity_info_t = flecsi::coloring::entity_info_t;
  
  // load the mesh
  auto filename_string = std::string( filename.data() );
  exodus_definition_t md( filename_string );

  // Create a communicator instance to get neighbor information.
  auto communicator = std::make_unique<flecsi::coloring::mpi_communicator_t>();
  auto comm_size = communicator->size();
  auto rank = communicator->rank();

  //----------------------------------------------------------------------------
  // Cell Coloring
  //----------------------------------------------------------------------------
    
  // Cells index coloring.
  flecsi::coloring::index_coloring_t cells;
  flecsi::coloring::coloring_info_t cell_color_info;
 
  // Create the dCRS representation for the distributed colorer.
  // This essentialy makes the graph of the dual mesh.
  auto dcrs = flecsi::coloring::make_dcrs(md);

  // Create a colorer instance to generate the primary coloring.
  auto colorer = std::make_unique<flecsi::coloring::parmetis_colorer_t>();

  // Create the primary coloring.
  cells.primary = colorer->color(dcrs);

  {
    clog_tag_guard(coloring);
    clog_container_one(info, "primary coloring", cells.primary, clog::space);
  } // guard

  //----------------------------------------------------------------------------
  // Cell Closure.  However many layers of ghost cells are needed are found 
  // here.
  //----------------------------------------------------------------------------
    
  // Compute the dependency closure of the primary cell coloring
  // through vertex intersections (specified by last argument "1").
  // To specify edge or face intersections, use 2 (edges) or 3 (faces).
  auto closure = 
    flecsi::topology::entity_neighbors<cell_dim, cell_dim, thru_dim>(
      md, cells.primary
    );

  {
    clog_tag_guard(coloring);
    clog_container_one(info, "closure", closure, clog::space);
  } // guard

  // Subtracting out the initial set leaves just the nearest
  // neighbors. This is similar to the image of the adjacency
  // graph of the initial indices.
  auto nearest_neighbors = 
    flecsi::utils::set_difference(closure, cells.primary);

  {
    clog_tag_guard(coloring);
    clog_container_one(
      info, "nearest neighbors", nearest_neighbors, clog::space
    );
  } // guard

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
      md, nearest_neighbors
    );

  {
    clog_tag_guard(coloring);
    clog_container_one(
      info, "nearest neighbor closure", nearest_neighbor_closure, clog::space
    );
  } // guard

  // Subtracting out the closure leaves just the
  // next nearest neighbors.
  auto next_nearest_neighbors =
    flecsi::utils::set_difference(nearest_neighbor_closure, closure);

  {
    clog_tag_guard(coloring);
    clog_container_one(
      info, "next nearest neighbor", next_nearest_neighbors, clog::space
    );
  } // guard

  // The union of the nearest and next-nearest neighbors gives us all
  // of the cells that might reference a vertex that we need.
  auto all_neighbors = 
    flecsi::utils::set_union(nearest_neighbors, next_nearest_neighbors);

  {
    clog_tag_guard(coloring);
    clog_container_one(info, "all neighbors", all_neighbors, clog::space);
  } // guard

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
      auto entry = entity_info_t(primary_indices_map[offset], rank, offset, i);
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

  {
    clog_tag_guard(coloring);
    clog_container_one(
      info, "exclusive cells ", cells.exclusive, clog::newline
    );
    clog_container_one(info, "shared cells ", cells.shared, clog::newline);
    clog_container_one(info, "ghost cells ", cells.ghost, clog::newline);
  } // guard
  
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
  
  flecsi::coloring::index_coloring_t vertices;
  flecsi::coloring::coloring_info_t vertex_color_info;

  color_entity<0>( 
    md, 
    communicator.get(), 
    closure, 
    remote_info_map, 
    shared_cells_map,
    closure_intersection_map,
    vertices, 
    vertex_color_info 
  );

  //----------------------------------------------------------------------------
  // Edge Closure
  //----------------------------------------------------------------------------
  
  flecsi::coloring::index_coloring_t edges;
  flecsi::coloring::coloring_info_t edge_color_info;

  color_entity<1>( 
    md, 
    communicator.get(), 
    closure, 
    remote_info_map, 
    shared_cells_map,
    closure_intersection_map,
    edges, 
    edge_color_info 
  );


  //----------------------------------------------------------------------------
  // Add the results to the context
  //----------------------------------------------------------------------------

  // Get the context instance.
  auto & context = flecsi::execution::context_t::instance();

  {
    clog_tag_guard(coloring);
    clog(info) << cell_color_info << std::endl << std::flush;
    clog(info) << vertex_color_info << std::endl << std::flush;
  } // gaurd

  // Gather the coloring info from all colors
  auto cell_coloring_info = 
    communicator->gather_coloring_info(cell_color_info);
  auto edge_coloring_info =
    communicator->gather_coloring_info(edge_color_info);
  auto vertex_coloring_info =
    communicator->gather_coloring_info(vertex_color_info);

  {
    clog_tag_guard(coloring);
    clog(info) << "vertex input coloring info color " 
               << rank << vertex_color_info << std::endl;
    for(auto ci: vertex_coloring_info) 
      clog(info) << "vertex coloring info color " << ci.first 
                 << ci.second << std::endl;
  }

  // Add colorings to the context.
  using index_spaces = mesh_t::index_spaces_t;

  context.add_coloring(index_spaces::vertices, vertices, vertex_coloring_info);
  context.add_coloring(index_spaces::edges, edges, edge_coloring_info);
  context.add_coloring(index_spaces::cells, cells, cell_coloring_info);
    
  //----------------------------------------------------------------------------
  // Add adjacency information
  //----------------------------------------------------------------------------
  
  // create the starting index
  size_t index = index_spaces::vertices_to_edges;

  // create a master list of all entities
  std::vector< std::vector<size_t> > entity_ids(num_dims+1);

  auto num_verts = 
    vertices.exclusive.size() + vertices.shared.size() + vertices.ghost.size();
  auto num_edges =
    edges.exclusive.size() + edges.shared.size() + edges.ghost.size();
  auto num_cells =
    cells.exclusive.size() + cells.shared.size() + cells.ghost.size();

  entity_ids[0].reserve( num_verts );
  entity_ids[1].reserve( num_edges );
  entity_ids[2].reserve( num_cells );

  for ( auto v : vertices.exclusive ) entity_ids[0].push_back(v.id);
  for ( auto v : vertices.shared    ) entity_ids[0].push_back(v.id);
  for ( auto v : vertices.ghost     ) entity_ids[0].push_back(v.id);

  for ( auto e : edges.exclusive ) entity_ids[1].push_back(e.id);
  for ( auto e : edges.shared    ) entity_ids[1].push_back(e.id);
  for ( auto e : edges.ghost     ) entity_ids[1].push_back(e.id);

  for ( auto c : cells.exclusive ) entity_ids[2].push_back(c.id);
  for ( auto c : cells.shared    ) entity_ids[2].push_back(c.id);
  for ( auto c : cells.ghost     ) entity_ids[2].push_back(c.id);

  // loop over each dimension and determine the adjacency sizes
  for ( int from_dim = 0; from_dim<=num_dims; ++from_dim ) {
   
    // the master list of all entity ids
    auto & from_ids = entity_ids[from_dim];

    // remove duplicate entries from the master lists
    std::sort( from_ids.begin(), from_ids.end() );
    from_ids.erase( 
      std::unique( from_ids.begin(), from_ids.end() ), 
      from_ids.end() 
    );

    for ( int to_dim = 0; to_dim<=num_dims; ++to_dim ) {

      // skip the case where both dimensions are the same
      if ( from_dim == to_dim ) continue;

      // populate the adjacency information
      flecsi::coloring::adjacency_info_t ai;
      ai.index_space = index++;
      ai.from_index_space = from_dim;
      ai.to_index_space = to_dim;
      ai.color_sizes.resize(comm_size);
  
      // loop over all cells and count the number of adjacencies
      size_t cnt = 0;
      for ( auto c : from_ids )
        cnt += md.entities(from_dim, to_dim, c).size();

      // gather the results
      ai.color_sizes = communicator->gather_sizes( cnt );
    
      // add the result to the context
      context.add_adjacency(ai);

    }
  }

  //----------------------------------------------------------------------------
  // output the result
  //----------------------------------------------------------------------------

  // some aliases
  using point_t = mesh_t::point_t; 
  using exodus_t = flecsi::io::exodus_base__<num_dims, real_t>;
  
  //------------------------------------
  // Open the file

  // figure out this ranks file name
  auto basename = utils::basename( filename_string );
  auto output_prefix = utils::remove_extension( basename );
  std::stringstream output_filename;
  output_filename << output_prefix;
  output_filename << "_rank";
  output_filename << std::setfill('0') << std::setw(6) << rank;
  output_filename << ".exo";

  // open the exodus file
  auto exoid = exodus_t::open( output_filename.str(), std::ios_base::out );

  //------------------------------------
  // Set exodus parameters

  // set exodus parameters
  auto num_nodes = md.num_entities( 0 );
  auto num_faces = num_dims==3 ? md.num_entities( num_dims-1 ) : 0;
  auto num_elems = 
    cells.exclusive.size() + cells.shared.size() + cells.ghost.size();;

  auto exo_params = exodus_t::make_params();
  exo_params.num_nodes = num_nodes;
  if ( num_dims == 3 ) {
    exo_params.num_face = num_faces;
    exo_params.num_face_blk = 1;
  }
  exo_params.num_elem = num_elems;
  exo_params.num_elem_blk = 
    !cells.exclusive.empty() + !cells.shared.empty() + !cells.ghost.empty();
  exo_params.num_node_sets = 
    !vertices.exclusive.empty() + !vertices.shared.empty() + 
    !vertices.ghost.empty();

  exodus_t::write_params(exoid, exo_params);

  //------------------------------------
  // Write the coordinates

  std::vector<real_t> vertex_coord( num_nodes * num_dims );
  
  for ( size_t i=0; i<num_nodes; ++i ) {
    const auto & vert = md.vertex<point_t>(i);
    for ( int d=0; d<num_dims; ++d ) 
      vertex_coord[ d*num_nodes + i ] = vert[d];
  }

  exodus_t::write_point_coords( exoid, vertex_coord );
  
  //------------------------------------
  // Write the faces
  
  if ( num_dims == 3 ) {
    exodus_t::template write_face_block<int>( 
      exoid, 1, "faces", md.entities(2,0)
    );
  }

  //------------------------------------
  // Write Exclusive Cells / vertices
  
  // the block id and side set counter
  int elem_blk_id = 0;
  int node_set_id = 0;

  // 3d wants cell faces, 2d wants cell vertices
  auto to_dim = (num_dims == 3) ? 2 : 0;
  const auto & cell_entities = md.entities(cell_dim, to_dim);

  // lambda function to convert to integer lists
  auto to_list = [&](const auto & list_in)
    -> std::vector<std::vector<size_t>>
  {
    std::vector<std::vector<size_t>> list_out;
    list_out.reserve( list_in.size() );
    for ( auto & e : list_in )
      list_out.emplace_back( cell_entities[e.id] );
    return list_out;
  };

  // write the cells
  exodus_t::template write_element_block<int>( 
    exoid, ++elem_blk_id, "exclusive cells", to_list(cells.exclusive) 
  );
  exodus_t::template write_element_block<int>( 
    exoid, ++elem_blk_id, "shared cells", to_list(cells.shared)
  );
  exodus_t::template write_element_block<int>( 
    exoid, ++elem_blk_id, "ghost cells", to_list(cells.ghost)
  );
    
  // lambda function to convert to integer lists
  auto to_vec = [](const auto & list_in)
    -> std::vector<size_t>
  {
    std::vector<size_t> list_out;
    list_out.reserve( list_in.size() );
    for ( auto & e : list_in )
      list_out.push_back( e.id );
    return list_out;
  };

  // write the vertices
  exodus_t::template write_node_set<int>(
    exoid, ++node_set_id, "exclusive vertices", to_vec(vertices.exclusive)
  );
  exodus_t::template write_node_set<int>( 
    exoid, ++node_set_id, "shared vertices", to_vec(vertices.shared)
  );
  exodus_t::template write_node_set<int>( 
    exoid, ++node_set_id, "ghost vertices", to_vec(vertices.ghost)
  );

      
  //------------------------------------
  // Close the file

  exodus_t::close( exoid );

  clog(info) << "Finished mesh partitioning." << std::endl;

} // somerhing

////////////////////////////////////////////////////////////////////////////////
/// \brief the main mesh initialization driver
////////////////////////////////////////////////////////////////////////////////
void initialize_mesh( 
  client_handle__<mesh_t, flecsi::dwd> mesh, 
  char_array_t filename
) {
  clog(info) << "INIT MESH TASK" << std::endl;

  // some constant expressions
  constexpr auto num_dims = mesh_t::num_dimensions;
  
  // alias some types
  using real_t = mesh_t::real_t;
  using exodus_definition_t = flecsi::io::exodus_definition__<num_dims, real_t>;
  using point_t = typename mesh_t::point_t;
  using vertex_t = typename mesh_t::vertex_t;
  using cell_t = typename mesh_t::cell_t;
  
  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // get the entity maps
  auto & vertex_map = context.index_map( mesh_t::index_spaces_t::vertices );
  auto & reverse_vertex_map = 
    context.reverse_index_map( mesh_t::index_spaces_t::vertices );

  auto & cell_map = context.index_map( mesh_t::index_spaces_t::cells );
  
  //----------------------------------------------------------------------------
  // Load the mesh
  //----------------------------------------------------------------------------
  
  auto filename_string = std::string( filename.data() );
  exodus_definition_t md( filename_string );

  //----------------------------------------------------------------------------
  // create the vertices
  //----------------------------------------------------------------------------

  std::vector< vertex_t * > vertices;
  vertices.reserve( vertex_map.size() );

  for(auto & vm: vertex_map) { 
    // get the point
    const auto & p = md.vertex<point_t>( vm.second );
    // now create it
    auto v = mesh.create_vertex( p );
    vertices.emplace_back(v);
  } // for vertices

  std::cout << "vertices: " << vertices.size() << std::endl;
  
  //----------------------------------------------------------------------------
  // create the cells
  //----------------------------------------------------------------------------

  for(auto & cm: cell_map) { 
    // get the list of vertices
    auto vs = 
      md.entities( cell_t::dimension, vertex_t::dimension, cm.second );
    // create a list of vertex pointers
    std::vector< vertex_t * > elem_vs( vs.size() );
    // transform the list of vertices to mesh ids
    std::transform(
      vs.begin(),
      vs.end(),
      elem_vs.begin(),
      [&](auto v) { return vertices[ reverse_vertex_map[v] ]; }
    );
    // create the cell
    auto c = mesh.create_cell( elem_vs ); 
  }

  std::cout << "cells: " << cell_map.size() << std::endl;
  
  //----------------------------------------------------------------------------
  // initialize the mesh
  //----------------------------------------------------------------------------
  
  mesh.init(); 

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
  enable_exceptions();

  //===========================================================================
  // Parse arguments
  //===========================================================================
  
  auto args = process_arguments( argc, argv );
  
  // process the simple ones
  if ( args.count("h") )
    return 0;
 
  // get the input file
  auto mesh_file_name = 
    args.count("m") ? args.at("m") : std::string();

  // override any inputs if need be
  if ( !mesh_file_name.empty() )
    std::cout << "Using mesh file \"" << mesh_file_name << "\"." 
              << std::endl;
  else
    raise_runtime_error( "No mesh file provided" );
  
  //===========================================================================
  // Partition mesh
  //===========================================================================
  
  clog(info) << "Partitioning mesh" << std::endl;
  
  // need to put the filename into a statically sized character array
  char_array_t filename;
  strcpy( filename.data(), mesh_file_name.c_str() );

  flecsi_execute_mpi_task(partition_mesh, filename);
  
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
  auto mesh_file_name = 
    args.count("m") ? args.at("m") : std::string();
  
  // need to put the filename into a statically sized character array
  char_array_t filename;
  strcpy( filename.data(), mesh_file_name.c_str() );
  
  //===========================================================================
  // Load the mesh
  //===========================================================================

  // get a mesh handle and call the initialization task
  auto mesh_handle = flecsi_get_client_handle(mesh_t, meshes, mesh0);
  auto f1 = flecsi_execute_task(initialize_mesh, single, mesh_handle, filename);
  f1.wait();

  return 0;
}


} // namespace
} // namespace
