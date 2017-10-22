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
#include <cinch/logging/cinchlog.h>
#include <flecsale/utils/trivial_string.h>
#include <flecsi/coloring/dcrs_utils.h>
#include <flecsi/coloring/mpi_communicator.h>
#include <flecsi/coloring/parmetis_colorer.h>
#include <flecsi/execution/execution.h>
#include <flecsi/io/exodus_definition.h>
#include <flecsi/topology/closure_utils.h>

// system includes
#include <iostream>

namespace apps {
namespace common {
  
////////////////////////////////////////////////////////////////////////////////
//! Build the list of exclusve, shared and ghost entities
////////////////////////////////////////////////////////////////////////////////
template < 
  typename COMM,
  typename CLOSURE_SET,
  typename ENTITY_MAP,
  typename INTERSECTION_MAP,
  typename CONNECTIVITY_MAP,
  typename INDEX_COLOR,
  typename COLOR_INFO
>
auto color_entity(
  COMM * communicator,
  const CLOSURE_SET & cells_plus_halo, 
  const ENTITY_MAP & remote_info_map,
  const ENTITY_MAP & shared_cells_map,
  const INTERSECTION_MAP & closure_intersection_map,
  const CONNECTIVITY_MAP & cell_to_entity_map,
  const CONNECTIVITY_MAP & entity_to_cell_map,
  INDEX_COLOR & entities, 
  COLOR_INFO & entity_color_info
) {
  
  // some type aliases
  using entity_info_t = flecsi::coloring::entity_info_t;

  // info about the mpi communicator
  auto comm_size = communicator->size();
  auto rank = communicator->rank();

  //----------------------------------------------------------------------------

  // Form the entity closure.  i.e. get a list of all entities that are
  // referenced by this rank's cells
  std::vector< size_t > referenced_entities;
  // guess an initial size
  referenced_entities.reserve( cells_plus_halo.size() );
  
  // now add all elements
  for(auto c: cells_plus_halo) {
    const auto & ents = cell_to_entity_map.at(c);
    referenced_entities.insert(
      referenced_entities.end(), ents.begin(), ents.end() );
  }

  // remove non-unique entries
  std::sort( referenced_entities.begin(), referenced_entities.end() );
  auto last = 
    std::unique( referenced_entities.begin(), referenced_entities.end() );
  referenced_entities.erase( last, referenced_entities.end() );
    

  //----------------------------------------------------------------------------
  
  // Assign entity ownership
  std::vector<std::set<size_t>> entity_requests(comm_size);
  std::set<entity_info_t> entity_info;

  size_t offset(0);
  for(auto i: referenced_entities) {

    // Get the set of cells that reference this entity.
    const auto & referencers = entity_to_cell_map.at(i);

    size_t min_rank(std::numeric_limits<size_t>::max());
    std::set<size_t> shared_entities;

    // Iterate the direct referencers to assign entity ownership.
    for(auto c: referencers) {

      // Check the remote info map to see if this cell is
      // off-color. If it is, compare it's rank for
      // the ownership logic below.
      auto it = remote_info_map.find(c);
      if(it != remote_info_map.end()) {
        min_rank = std::min(min_rank, it->second.rank);
        shared_entities.insert(it->second.rank);
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
  
  entity_color_info.exclusive = entities.exclusive.size();
  entity_color_info.shared = entities.shared.size();
  entity_color_info.ghost = entities.ghost.size();
}


////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 2D version
////////////////////////////////////////////////////////////////////////////////
template<
  typename MESH_DEFINITION,
  typename MESH_TYPE,
  bool Enabled = ( std::decay_t<MESH_DEFINITION>::dimension() == 2 ),
  typename = std::enable_if_t< Enabled >
>
void create_cells( MESH_DEFINITION && mesh_def, MESH_TYPE && mesh )
{
  
  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------
 
  using mesh_t = typename std::decay_t< MESH_TYPE >;

  // some constant expressions
  constexpr auto num_dims = mesh_t::num_dimensions;
  
  // Alias the index spaces type
  using index_spaces = typename mesh_t::index_spaces_t;

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
  typename MESH_DEFINITION,
  typename MESH_TYPE,
  bool Enabled = ( std::decay_t<MESH_DEFINITION>::dimension() == 3 ),
  typename std::enable_if_t< Enabled >* = nullptr
>
void create_cells( MESH_DEFINITION && mesh_def, MESH_TYPE && mesh )
{
  
  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------
  
  using mesh_t = typename std::decay_t< MESH_TYPE >;
  
  // some constant expressions
  constexpr auto num_dims = mesh_t::num_dimensions;
  
  // Alias the index spaces type
  using index_spaces = typename mesh_t::index_spaces_t;

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
/// \brief Helper to make corners
////////////////////////////////////////////////////////////////////////////////
template<
  typename MESH_DEFINITION
>
auto make_corners( const MESH_DEFINITION & mesh_def )
{
  // not sure why this needs to be decayed
  using mesh_definition_t = std::decay_t<MESH_DEFINITION>;
  using connectivity_t = typename mesh_definition_t::connectivity_t;

  // get the number of dimensions
  constexpr auto num_dims = mesh_definition_t::dimension();

  // get the connectvitiy from the cells to the vertices.  there is a 
  // corner per cell, per vertex
  const auto & cells_to_vertices = mesh_def.entities(num_dims, 0);
  auto num_cells = cells_to_vertices.size();

  // count corners
  size_t num_corners = 0;
  for ( const auto & verts : cells_to_vertices )
    num_corners += verts.size();
  
  // create storage
  connectivity_t cells_to_corners;
  connectivity_t corners_to_cells;

  cells_to_corners.resize( num_cells );
  corners_to_cells.resize( num_corners );
  
  // now make them
  size_t corner_id = 0;
  size_t cell_id = 0;

  for ( const auto & verts : cells_to_vertices ) {
    std::vector< size_t > corner_ids;
    corner_ids.reserve( verts.size() );
    for ( auto v : verts ) {
      corner_ids.emplace_back( corner_id );
      corners_to_cells[corner_id].emplace_back( cell_id );
      ++corner_id;
    }
    cells_to_corners[cell_id] = corner_ids;
    ++cell_id;
  }

  return std::make_pair( cells_to_corners, corners_to_cells );

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper to make wedges
////////////////////////////////////////////////////////////////////////////////
template<
  typename MESH_DEFINITION
>
auto make_wedges( const MESH_DEFINITION & mesh_def )
{

  // not sure why this needs to be decayed
  using mesh_definition_t = std::decay_t<MESH_DEFINITION>;
  using connectivity_t = typename mesh_definition_t::connectivity_t;

  // get the number of dimensions
  constexpr auto num_dims = mesh_definition_t::dimension();

  // get the connectvitiy from the cells to the vertices.  there is a 
  // wedge per cell, per vertex
  const auto & cells_to_faces = mesh_def.entities(num_dims, num_dims-1);
  // const auto & faces_to_edges = mesh_def.entities(num_dims-1, num_dims-2);
  auto num_cells = mesh_def.num_entities(num_dims);

  // count wedges
  size_t num_wedges = 0;
  for ( const auto & faces : cells_to_faces )
    num_wedges += 2*faces.size();
  
  // create storage
  connectivity_t cells_to_wedges;
  connectivity_t wedges_to_cells;

  cells_to_wedges.resize( num_cells );
  wedges_to_cells.resize( num_wedges );
  
  // now make them
  size_t wedge_id = 0;
  size_t cell_id = 0;

  for ( const auto & faces : cells_to_faces ) {
    std::vector< size_t > wedge_ids;
    wedge_ids.reserve( (num_dims+1) * faces.size() ); // an estimate for size
    for ( auto f : faces ) {
      //for ( auto e : faces_to_edges.at(f) ) {
        // two wedges for each edge ( one attached to each vertex )
        for ( int i=0; i<2; ++i ) {
          wedge_ids.emplace_back( wedge_id );
          wedges_to_cells[wedge_id].emplace_back( cell_id );
          ++wedge_id;
        }
      //}
    }
    cells_to_wedges[cell_id] = wedge_ids;
    ++cell_id;
  }

  return std::make_pair( cells_to_wedges, wedges_to_cells );

}

} // namespace
} // namespace
