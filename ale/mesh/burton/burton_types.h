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
 * \file burton_types.h
 * \authors bergen
 * \date Initial file creation: Sep 02, 2015
 ******************************************************************************/

#pragma once

//! user incldues
#include "flecsi/mesh/mesh_topology.h"
#include "../../mesh/burton/burton_mesh_traits.h"
#include "../../mesh/burton/burton_entity_types.h"

namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_mesh_types_t burton_types.h
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_mesh_types_t {

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using mesh_traits_t = burton_mesh_traits_t<N>;

  //! The dimension of the burton mesh picked up from burton_mesh_traits_t.
  static constexpr size_t dimension = mesh_traits_t::dimension;

  //! The number of domains in burton mesh picked up from burton_mesh_traits_t.
  static constexpr size_t num_domains = mesh_traits_t::num_domains;

  //============================================================================
  // Define basic types.
  //============================================================================

  //! Type for burton mesh vertices.
  using vertex_t = burton_vertex_t<N>;

  //! Type for burton mesh edges.
  using edge_t = burton_edge_t<N>;

  // Cell types
  //! Type for burton mesh cells.
  using cell_t = burton_cell_t<N>;

  //! Types for burton mesh cells.
  using triangle_cell_t = burton_triangle_cell_t<N>;
  using quadrilateral_cell_t = burton_quadrilateral_cell_t<N>;
  using polygonal_cell_t = burton_polygonal_cell_t<N>;

  //! Type for burton mesh corners.
  using corner_t = burton_corner_t<N>;

  //! Type for burton mesh wedges.
  using wedge_t = burton_wedge_t<N>;

  //============================================================================
  // Specify mesh parameterizations.
  //============================================================================

  //! Definitions of burton mesh entities and their domain.
  // clang-format off
  using entity_types =
      std::tuple<
        std::pair<flecsi::domain_<0>, vertex_t>,
        std::pair<flecsi::domain_<0>, edge_t>,
        std::pair<flecsi::domain_<0>, cell_t>,
        std::pair<flecsi::domain_<1>, wedge_t>,
        std::pair<flecsi::domain_<1>, corner_t>
      >;

  //! Connectivities are adjacencies of entities within a single domain.
  using connectivities =
    std::tuple<
      std::tuple<flecsi::domain_<0>, vertex_t, edge_t>,
      std::tuple<flecsi::domain_<0>, vertex_t, cell_t>,
      std::tuple<flecsi::domain_<0>, edge_t, vertex_t>,
      std::tuple<flecsi::domain_<0>, edge_t, cell_t>,
      std::tuple<flecsi::domain_<0>, cell_t, vertex_t>,
      std::tuple<flecsi::domain_<0>, cell_t, edge_t>
      >;

  //! Bindings are adjacencies of entities across two domains.
  using bindings =
      std::tuple<
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>, cell_t, corner_t>,
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>, edge_t, corner_t>,
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>, vertex_t, corner_t>,
        std::tuple<flecsi::domain_<1>, flecsi::domain_<0>, corner_t, cell_t>,
        std::tuple<flecsi::domain_<1>, flecsi::domain_<0>, corner_t, edge_t>,
        std::tuple<flecsi::domain_<1>, flecsi::domain_<0>, corner_t, vertex_t>
      >;
  // clang-format on

}; // struct burton_mesh_types_t

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
