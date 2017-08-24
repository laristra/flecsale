/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file 
/// \brief Associate the mesh entities with flecsi.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user incldues
#include "flecsale/geom/shapes/geometric_shapes.h"

#include "flecsale/mesh/burton/burton_vertex.h"
#include "flecsale/mesh/burton/burton_corner.h"
#include "flecsale/mesh/burton/burton_element.h"
#include "flecsale/mesh/burton/burton_wedge.h"

#include "flecsi/topology/mesh.h"
#include "flecsale/mesh/burton/burton_config.h"

namespace flecsale {
namespace mesh {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
//! \tparam N  The number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N, bool Include_Extras = false >
struct burton_types_t {};

////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
//! \remark This is the two-dimensional version.
////////////////////////////////////////////////////////////////////////////////
struct burton_2d_types_base
{

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using config_t = burton_config_t<2>;

  //! The dimension of the burton mesh picked up from burton_config_t.
  static constexpr size_t num_dimensions = config_t::num_dimensions;

  //! The number of domains in burton mesh picked up from burton_config_t.
  static constexpr size_t num_domains = config_t::num_domains;

  //! the flecsi mesh topology type
  using mesh_topology_base_t = typename config_t::mesh_topology_base_t;

  //! the base type for the entities
  using mesh_entity_base_t = typename config_t::mesh_entity_base_t;

  //! the id type
  using id_t = flecsi::utils::id_t;

  //============================================================================
  //! \brief The burton mesh index spaces.
  //============================================================================
  struct index_spaces_t {

    //! The individual enumeration of the index spaces
    enum index_spaces : size_t {
      // the main index spaces
      vertices,
      edges,
      faces = edges,
      cells,
      corners,
      wedges,
      // index spaces for connectivity
      vertices_to_edges,
      vertices_to_faces = vertices_to_edges,
      vertices_to_cells,
      edges_to_vertices,
      edges_to_cells,
      faces_to_vertices = edges_to_vertices,
      faces_to_cells = edges_to_cells,
      cells_to_vertices,
      cells_to_edges,
      cells_to_faces = cells_to_edges,
      // index spaces that are not used
      edges_to_faces = 7777,
      faces_to_edges = 7777
    };

    //! Maps an entity dimension to an index space id
    static constexpr size_t entity_map[2][3] = {
      // domain
      vertices,
      edges,
      cells,
      // domain
      corners,
      wedges,
      7777
    };


    //! Maps dimension-to-dimension connectivity to an index space id
    static constexpr size_t connectivity_map[3][3] = {
      // row
      7777,
      vertices_to_edges,
      vertices_to_cells,
      // row
      edges_to_vertices,
      7777,
      edges_to_cells,
      // row
      cells_to_vertices,
      cells_to_edges,
      7777
    };


  };
 

  //============================================================================
  // Define basic types.
  //============================================================================

  //! Type for burton mesh vertices.
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! Type for burton mesh edges.
  using edge_t = burton_edge_t<num_dimensions>;

  //! Type for burton mesh faces.
  using face_t = burton_face_t<num_dimensions>;

  //! Type for burton mesh cells.
  using cell_t = burton_cell_t<num_dimensions>;

  //! Type for burton mesh corners.
  using corner_t = burton_corner_t<num_dimensions>;

  //! Type for burton mesh wedges.
  using wedge_t = burton_wedge_t<num_dimensions>;


  //============================================================================
  //! \brief depending upon the dimension/number of verices, create different 
  //!   types of entities
  //! \tparam M The domain index.
  //! \tparam D The dimensional index.
  //============================================================================
  template<size_t M, size_t D>
  static constexpr 
  mesh_entity_base_t *
  create_entity(mesh_topology_base_t* mesh, size_t num_vertices, const id_t & id) 
  {
    switch(M){
      //---------- Primal Mesh ----------//
    case 0:
      switch(D) {
      case 1:
        return mesh->template make<edge_t>(id);
      default:
        raise_logic_error("invalid topological dimension");
      }
      //---------- Dual Mesh ----------//
    case 1:
      switch(D) {
      case 0:
        return mesh->template make<corner_t>();
      case 1:
        return mesh->template make<wedge_t>();
      default:
        raise_logic_error("invalid topological dimension");
      }
      //---------- Error ----------//
    default:
      raise_logic_error("invalid domain");
    }
    // should never get here
    return nullptr;
  }
  
};

////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
//! \remark This is the two-dimensional version.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_types_t<2, false> : public burton_2d_types_base
{

  //! Definitions of burton mesh entities and their domain.
  flecsi_register_entity_types(
    flecsi_entity_type( index_spaces_t::vertices, 0, vertex_t ),
    flecsi_entity_type( index_spaces_t::edges, 0, edge_t ),
    flecsi_entity_type( index_spaces_t::cells, 0, cell_t )
  );


  //! Connectivities are adjacencies of entities within a single domain.
  flecsi_register_connectivities(
    flecsi_connectivity( index_spaces_t::vertices_to_edges, 0, vertex_t, edge_t ),
    flecsi_connectivity( index_spaces_t::vertices_to_cells, 0, vertex_t, cell_t ),
    flecsi_connectivity( index_spaces_t::edges_to_vertices, 0, edge_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::edges_to_cells,    0, edge_t, cell_t ),
    flecsi_connectivity( index_spaces_t::cells_to_vertices, 0, cell_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::cells_to_edges,    0, cell_t, edge_t )
  );

  //! Bindings are adjacencies of entities across two domains.
  flecsi_register_bindings();

};

////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
//! \remark This is the two-dimensional version.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_types_t<2, true> : public burton_2d_types_base
{

  //! Definitions of burton mesh entities and their domain.
  flecsi_register_entity_types(
    flecsi_entity_type( index_spaces_t::vertices, 0, vertex_t ),
    flecsi_entity_type( index_spaces_t::edges, 0, edge_t ),
    flecsi_entity_type( index_spaces_t::cells, 0, cell_t ),
    flecsi_entity_type( index_spaces_t::corners, 1, corner_t ),
    flecsi_entity_type( index_spaces_t::wedges, 1, wedge_t )
  );


  //! Connectivities are adjacencies of entities within a single domain.
  flecsi_register_connectivities(
    flecsi_connectivity( index_spaces_t::vertices_to_edges, 0, vertex_t, edge_t ),
    flecsi_connectivity( index_spaces_t::vertices_to_cells, 0, vertex_t, cell_t ),
    flecsi_connectivity( index_spaces_t::edges_to_vertices, 0, edge_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::edges_to_cells,    0, edge_t, cell_t ),
    flecsi_connectivity( index_spaces_t::cells_to_vertices, 0, cell_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::cells_to_edges,    0, cell_t, edge_t )
  );

  //! Bindings are adjacencies of entities across two domains.
  flecsi_register_bindings();
#if 0
     // corners
     std::tuple<index_space_<12>, domain_<0>, domain_<1>, cell_t,   corner_t>,
     std::tuple<index_space_<13>, domain_<0>, domain_<1>, edge_t,   corner_t>,
     std::tuple<index_space_<14>, domain_<0>, domain_<1>, vertex_t, corner_t>,
     std::tuple<index_space_<15>, domain_<1>, domain_<0>, corner_t,   cell_t>,
     std::tuple<index_space_<16>, domain_<1>, domain_<0>, corner_t,   edge_t>,
     std::tuple<index_space_<17>, domain_<1>, domain_<0>, corner_t, vertex_t>,
     // wedges
     std::tuple<index_space_<18>, domain_<0>, domain_<1>, cell_t,    wedge_t>,
     std::tuple<index_space_<19>, domain_<0>, domain_<1>, edge_t,    wedge_t>,
     std::tuple<index_space_<20>, domain_<0>, domain_<1>, vertex_t,  wedge_t>,
     std::tuple<index_space_<21>, domain_<1>, domain_<0>, wedge_t,    cell_t>,
     std::tuple<index_space_<22>, domain_<1>, domain_<0>, wedge_t,    edge_t>,
     std::tuple<index_space_<23>, domain_<1>, domain_<0>, wedge_t,  vertex_t>,
     // corner <-> wedges
     std::tuple<index_space_<24>, domain_<1>, domain_<1>,  wedge_t,  corner_t>,
     std::tuple<index_space_<25>, domain_<1>, domain_<1>,  corner_t, wedge_t>
  );
#endif

};


////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
struct burton_3d_types_base
{

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using config_t = burton_config_t<3>;

  //! The dimension of the burton mesh picked up from burton_config_t.
  static constexpr size_t num_dimensions = config_t::num_dimensions;

  //! The number of domains in burton mesh picked up from burton_config_t.
  static constexpr size_t num_domains = config_t::num_domains;

  //! the shape type
  using shape_t = config_t::shape_t;

  //! the flecsi mesh topology type
  using mesh_topology_base_t = typename config_t::mesh_topology_base_t;
  
  //! the base type for the entities
  using mesh_entity_base_t = typename config_t::mesh_entity_base_t;

  //! the id type
  using id_t = flecsi::utils::id_t;

  //============================================================================
  //! \brief The burton mesh index spaces.
  //============================================================================
  struct index_spaces_t {

    //! The individual enumeration of the index spaces
    enum index_spaces : size_t {
      // the main index spaces
      vertices,
      edges,
      faces,
      cells,
      corners,
      wedges,
      // index spaces for connectivity
      vertices_to_edges,
      vertices_to_faces,
      vertices_to_cells,
      edges_to_vertices,
      edges_to_faces,
      edges_to_cells,
      faces_to_vertices,
      faces_to_edges,
      faces_to_cells,
      cells_to_vertices,
      cells_to_edges,
      cells_to_faces,
    };

    //! Maps an entity dimension to an index space id
    static constexpr size_t entity_map[2][4] = {
      // domain
      vertices,
      edges,
      faces,
      cells,
      // domain
      corners,
      wedges,
      7777,
      7777
    };


    //! Maps dimension-to-dimension connectivity to an index space id
    static constexpr size_t connectivity_map[4][4] = {
      // row
      7777,
      vertices_to_edges,
      vertices_to_faces,
      vertices_to_cells,
      // row
      edges_to_vertices,
      7777,
      edges_to_faces,
      edges_to_cells,
      // row
      faces_to_vertices,
      faces_to_edges,
      7777,
      faces_to_cells,
      // row
      cells_to_vertices,
      cells_to_edges,
      cells_to_faces,
      7777
    };


  };
 

  //============================================================================
  // Define basic types.
  //============================================================================

  //! Type for burton mesh vertices.
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! Type for burton mesh edges.
  using edge_t = burton_edge_t<num_dimensions>;

  //! Type for burton mesh cells.
  using face_t = burton_face_t<num_dimensions>;

  //! Type for burton mesh cells.
  using cell_t = burton_cell_t<num_dimensions>;

  //! Type for burton mesh corners.
  using corner_t = burton_corner_t<num_dimensions>;

  //! Type for burton mesh wedges.
  using wedge_t = burton_wedge_t<num_dimensions>;

  //============================================================================
  //! \brief depending upon the dimension/number of verices, create different 
  //!   types of face entities
  //============================================================================  
  static
  mesh_entity_base_t *
  create_face(mesh_topology_base_t* mesh, size_t num_vertices, const id_t & id)
  {
    auto face_type = shape_t::polygon;
    switch(num_vertices) {
    case (1,2):
      raise_runtime_error( "can't have <3 vertices" );
    case (3):
      face_type = shape_t::triangle;
    case (4):
      face_type = shape_t::quadrilateral;
    }
    return mesh->template make<face_t>(id, face_type);
  }


  //============================================================================
  //! \brief depending upon the dimension/number of verices, create different 
  //!   types of entities
  //! \tparam M The domain index.
  //! \tparam D The dimensional index.
  //============================================================================
  template<size_t M, size_t D>
  static constexpr 
  mesh_entity_base_t *
  create_entity(mesh_topology_base_t* mesh, size_t num_vertices, const id_t & id) 
  {
    switch(M){
      //---------- Primal Mesh ----------//
    case 0:
      switch(D) {
      case 1:
        return mesh->template make<edge_t>(id);
      case 2:
        return create_face( mesh, num_vertices, id );
      default:
        raise_logic_error("invalid topological dimensions");
      }
#if 0
      //---------- Dual Mesh ----------//
    case 1:
      switch(D) {
      case 0:
        return mesh->template make<corner_t>();
      case 1:
        return mesh->template make<wedge_t>();
      default:
        raise_logic_error("invalid topological dimension");
      }
      //---------- Error ----------//
#endif
    default:
      raise_logic_error("invalid domain");
    }
    // should never get here
    return nullptr;
  }

};

////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_types_t<3, false> : public burton_3d_types_base
{
  
  //! Definitions of burton mesh entities and their domain.
  flecsi_register_entity_types(
    flecsi_entity_type( index_spaces_t::vertices, 0, vertex_t ),
    flecsi_entity_type( index_spaces_t::edges, 0, edge_t ),
    flecsi_entity_type( index_spaces_t::faces, 0, face_t ),
    flecsi_entity_type( index_spaces_t::cells, 0, cell_t )
    //flecsi_entity_type( attributes::wedges, 1, wedge_t ),
    //flecsi_entity_type( attributes::corners, 1, corner_t )
  );


  //! Connectivities are adjacencies of entities within a single domain.
  flecsi_register_connectivities(
    flecsi_connectivity( index_spaces_t::vertices_to_edges, 0, vertex_t, edge_t ),
    flecsi_connectivity( index_spaces_t::vertices_to_faces, 0, vertex_t, face_t ),
    flecsi_connectivity( index_spaces_t::vertices_to_cells, 0, vertex_t, cell_t ),
    flecsi_connectivity( index_spaces_t::edges_to_vertices, 0, edge_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::edges_to_faces,    0, edge_t,   face_t ),
    flecsi_connectivity( index_spaces_t::edges_to_cells,    0, edge_t,   cell_t ),
    flecsi_connectivity( index_spaces_t::faces_to_vertices, 0, face_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::faces_to_edges,    0, face_t,   edge_t ),
    flecsi_connectivity( index_spaces_t::faces_to_cells,    0, face_t,   cell_t ),
    flecsi_connectivity( index_spaces_t::cells_to_vertices, 0, cell_t, vertex_t ),
    flecsi_connectivity( index_spaces_t::cells_to_faces,    0, cell_t,   face_t ),
    flecsi_connectivity( index_spaces_t::cells_to_edges,    0, cell_t,   edge_t )
  );

  //! Bindings are adjacencies of entities across two domains.
  flecsi_register_bindings();
#if 0
  using bindings =
      std::tuple<
        // corners
        std::tuple<index_space_<18>, domain_<0>, domain_<1>,   cell_t, corner_t>,
        std::tuple<index_space_<19>, domain_<0>, domain_<1>,   face_t, corner_t>,
        std::tuple<index_space_<20>, domain_<0>, domain_<1>,   edge_t, corner_t>,
        std::tuple<index_space_<21>, domain_<0>, domain_<1>, vertex_t, corner_t>,
        std::tuple<index_space_<22>, domain_<1>, domain_<0>, corner_t,   cell_t>,
        std::tuple<index_space_<23>, domain_<1>, domain_<0>, corner_t,   face_t>,
        std::tuple<index_space_<24>, domain_<1>, domain_<0>, corner_t,   edge_t>,
        std::tuple<index_space_<25>, domain_<1>, domain_<0>, corner_t, vertex_t>,
        // wedges
        std::tuple<index_space_<26>, domain_<0>, domain_<1>,   cell_t,  wedge_t>,
        std::tuple<index_space_<27>, domain_<0>, domain_<1>,   face_t,  wedge_t>,
        std::tuple<index_space_<28>, domain_<0>, domain_<1>,   edge_t,  wedge_t>,
        std::tuple<index_space_<29>, domain_<0>, domain_<1>, vertex_t,  wedge_t>,
        std::tuple<index_space_<30>, domain_<1>, domain_<0>,  wedge_t,   cell_t>,
        std::tuple<index_space_<31>, domain_<1>, domain_<0>,  wedge_t,   face_t>,
        std::tuple<index_space_<32>, domain_<1>, domain_<0>,  wedge_t,   edge_t>,
        std::tuple<index_space_<33>, domain_<1>, domain_<0>,  wedge_t, vertex_t>,
        // corner <-> wedges
        std::tuple<index_space_<34>, domain_<1>, domain_<1>,  wedge_t,  corner_t>,
        std::tuple<index_space_<35>, domain_<1>, domain_<1>,  corner_t, wedge_t>
      >;
#endif

};

} // namespace burton
} // namespace mesh
} // namespace flecsale
