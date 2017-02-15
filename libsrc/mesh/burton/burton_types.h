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
#include "mesh/burton/burton_vertex.h"
#include "mesh/burton/burton_corner.h"
#include "mesh/burton/burton_element.h"
#include "mesh/burton/burton_wedge.h"

#include "mesh/burton/burton_triangle.h"
#include "mesh/burton/burton_quadrilateral.h"
#include "mesh/burton/burton_polygon.h"
#include "mesh/burton/burton_hexahedron.h"
#include "mesh/burton/burton_tetrahedron.h"
#include "mesh/burton/burton_polyhedron.h"

#include "flecsi/topology/mesh_topology.h"
#include "mesh/burton/burton_config.h"

namespace ale {
namespace mesh {



////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
//! \tparam N  The number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_types_t {};

////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
//! \remark This is the two-dimensional version.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_types_t<2> {

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using config_t = burton_config_t<2>;

  //! The dimension of the burton mesh picked up from burton_config_t.
  static constexpr size_t num_dimensions = config_t::num_dimensions;

  //! The number of domains in burton mesh picked up from burton_config_t.
  static constexpr size_t num_domains = config_t::num_domains;

  //! the base type for the mesh topology
  using mesh_topology_base_t = flecsi::topology::mesh_topology_base_t;

  // the base type for the entities
  using mesh_entity_base_t = flecsi::topology::mesh_entity_base_t<num_domains>;

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
  // Specify mesh parameterizations.
  //============================================================================

  //! Convenience type
  template<size_t D>
  using domain_ = flecsi::topology::domain_<D>;

  //! Definitions of burton mesh entities and their domain.
  using entity_types =
      std::tuple<
        std::pair<domain_<0>, vertex_t>,
        std::pair<domain_<0>, edge_t>,
        std::pair<domain_<0>, cell_t>,
        std::pair<domain_<1>, wedge_t>,
        std::pair<domain_<1>, corner_t>
      >;

  //! Connectivities are adjacencies of entities within a single domain.
  using connectivities =
    std::tuple<
      std::tuple<domain_<0>, vertex_t, edge_t>,
      std::tuple<domain_<0>, vertex_t, cell_t>,
      std::tuple<domain_<0>, edge_t, vertex_t>,
      // edges->edges makes sure edges(faces) works in 2d
      // std::tuple<domain_<0>, edge_t, edge_t>,
      std::tuple<domain_<0>, edge_t, cell_t>,
      std::tuple<domain_<0>, cell_t, vertex_t>,
      std::tuple<domain_<0>, cell_t, edge_t>
      >;

  //! Bindings are adjacencies of entities across two domains.
  using bindings =
      std::tuple<
        // corners
        std::tuple<domain_<0>, domain_<1>, cell_t,   corner_t>,
        std::tuple<domain_<0>, domain_<1>, edge_t,   corner_t>,
        std::tuple<domain_<0>, domain_<1>, vertex_t, corner_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t,   cell_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t,   edge_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t, vertex_t>,
        // wedges
        std::tuple<domain_<0>, domain_<1>, cell_t,    wedge_t>,
        std::tuple<domain_<0>, domain_<1>, edge_t,    wedge_t>,
        std::tuple<domain_<0>, domain_<1>, vertex_t,  wedge_t>,
        std::tuple<domain_<1>, domain_<0>, wedge_t,    cell_t>,
        std::tuple<domain_<1>, domain_<0>, wedge_t,    edge_t>,
        std::tuple<domain_<1>, domain_<0>, wedge_t,  vertex_t>,
        // corner <-> wedges
        std::tuple<domain_<1>, domain_<1>,  wedge_t,  corner_t>,
        std::tuple<domain_<1>, domain_<1>,  corner_t, wedge_t>
      >;


  //============================================================================
  //! \brief depending upon the dimension/number of verices, create different 
  //!   types of entities
  //! \tparam M The domain index.
  //! \tparam D The dimensional index.
  //============================================================================
  template<size_t M, size_t D>
  static constexpr 
  mesh_entity_base_t *
  create_entity(mesh_topology_base_t* mesh, size_t num_vertices) 
  {
    switch(M){
      //---------- Primal Mesh ----------//
    case 0:
      switch(D) {
      case 1:
        return mesh->template make<edge_t>(*mesh);
      default:
        raise_logic_error("invalid topological dimension");
      }
      //---------- Dual Mesh ----------//
    case 1:
      switch(D) {
      case 0:
        return mesh->template make<corner_t>(*mesh);
      case 1:
        return mesh->template make<wedge_t>(*mesh);
      default:
        raise_logic_error("invalid topological dimension");
      }
      //---------- Error ----------//
    default:
      raise_logic_error("invalid domain");
    }
  }
  
};


////////////////////////////////////////////////////////////////////////////////
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_types_t<3> {

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using config_t = burton_config_t<3>;

  //! The dimension of the burton mesh picked up from burton_config_t.
  static constexpr size_t num_dimensions = config_t::num_dimensions;

  //! The number of domains in burton mesh picked up from burton_config_t.
  static constexpr size_t num_domains = config_t::num_domains;

  //! the base type for the mesh topology
  using mesh_topology_base_t = flecsi::topology::mesh_topology_base_t;
  
  // the base type for the entities
  using mesh_entity_base_t = flecsi::topology::mesh_entity_base_t<num_domains>;
  
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
  // Specify mesh parameterizations.
  //============================================================================

  //! Convenience type
  template<size_t D>
  using domain_ = flecsi::topology::domain_<D>;

  //! Definitions of burton mesh entities and their domain.
  using entity_types =
      std::tuple<
        std::pair<domain_<0>, vertex_t>,
        std::pair<domain_<0>, edge_t>,
        std::pair<domain_<0>, face_t>,
        std::pair<domain_<0>, cell_t>,
        std::pair<domain_<1>, wedge_t>,
        std::pair<domain_<1>, corner_t>
      >;

  //! Connectivities are adjacencies of entities within a single domain.
  using connectivities =
    std::tuple<
      std::tuple<domain_<0>, vertex_t, edge_t>,
      std::tuple<domain_<0>, vertex_t, face_t>,
      std::tuple<domain_<0>, vertex_t, cell_t>,
      std::tuple<domain_<0>, edge_t, vertex_t>,
      std::tuple<domain_<0>, edge_t, face_t>,
      std::tuple<domain_<0>, edge_t, cell_t>,
      std::tuple<domain_<0>, face_t, vertex_t>,
      std::tuple<domain_<0>, face_t, edge_t>,
      std::tuple<domain_<0>, face_t, cell_t>,
      std::tuple<domain_<0>, cell_t, vertex_t>,
      std::tuple<domain_<0>, cell_t, face_t>,
      std::tuple<domain_<0>, cell_t, edge_t>
      >;

  //! Bindings are adjacencies of entities across two domains.
  using bindings =
      std::tuple<
        // corners
        std::tuple<domain_<0>, domain_<1>,   cell_t, corner_t>,
        std::tuple<domain_<0>, domain_<1>,   face_t, corner_t>,
        std::tuple<domain_<0>, domain_<1>,   edge_t, corner_t>,
        std::tuple<domain_<0>, domain_<1>, vertex_t, corner_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t,   cell_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t,   face_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t,   edge_t>,
        std::tuple<domain_<1>, domain_<0>, corner_t, vertex_t>,
        // wedges
        std::tuple<domain_<0>, domain_<1>,   cell_t,  wedge_t>,
        std::tuple<domain_<0>, domain_<1>,   face_t,  wedge_t>,
        std::tuple<domain_<0>, domain_<1>,   edge_t,  wedge_t>,
        std::tuple<domain_<0>, domain_<1>, vertex_t,  wedge_t>,
        std::tuple<domain_<1>, domain_<0>,  wedge_t,   cell_t>,
        std::tuple<domain_<1>, domain_<0>,  wedge_t,   face_t>,
        std::tuple<domain_<1>, domain_<0>,  wedge_t,   edge_t>,
        std::tuple<domain_<1>, domain_<0>,  wedge_t, vertex_t>,
        // corner <-> wedges
        std::tuple<domain_<1>, domain_<1>,  wedge_t,  corner_t>,
        std::tuple<domain_<1>, domain_<1>,  corner_t, wedge_t>
      >;
 
  //============================================================================
  //! \brief depending upon the dimension/number of verices, create different 
  //!   types of face entities
  //============================================================================  
  static
  mesh_entity_base_t *
  create_face(mesh_topology_base_t* mesh, size_t num_vertices) {
    switch(num_vertices) {
    case (1,2):
      raise_runtime_error( "can't have <3 vertices" );
    case (3):
      return mesh->template make< burton_triangle_t<num_dimensions> >(*mesh);
    case (4):
      return mesh->template make< burton_quadrilateral_t<num_dimensions> >(*mesh);
      break;
    default:
      return mesh->template make< burton_polygon_t<num_dimensions> >(*mesh);
      break;
    }
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
  create_entity(mesh_topology_base_t* mesh, size_t num_vertices)
  {   
    switch(M){
      //---------- Primal Mesh ----------//
    case 0:
      switch(D) {
      case 1:
        return mesh->template make<edge_t>(*mesh);
      case 2:
        return create_face(mesh, num_vertices);
      default:
        raise_logic_error("invalid topological dimensions");
      }
      //---------- Dual Mesh ----------//
    case 1:
      switch(D) {
      case 0:
        return mesh->template make<corner_t>(*mesh);
      case 1:
        return mesh->template make<wedge_t>(*mesh);
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

} // namespace mesh
} // namespace ale
