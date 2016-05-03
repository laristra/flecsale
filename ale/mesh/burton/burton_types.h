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
#include "ale/mesh/burton/burton_vertex.h"
#include "ale/mesh/burton/burton_corner.h"
#include "ale/mesh/burton/burton_element.h"
#include "ale/mesh/burton/burton_wedge.h"

#include "flecsi/mesh/mesh_topology.h"
#include "ale/mesh/burton/burton_mesh_traits.h"

namespace ale {
namespace mesh {



////////////////////////////////////////////////////////////////////////////////
//! \class burton_mesh_types_t burton_types.h
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_mesh_types_t {};

////////////////////////////////////////////////////////////////////////////////
//! \class burton_mesh_types_t burton_types.h
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_mesh_types_t<2> {

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using mesh_traits_t = burton_mesh_traits_t<2>;

  //! The dimension of the burton mesh picked up from burton_mesh_traits_t.
  static constexpr size_t dimension = mesh_traits_t::dimension;

  //! The number of domains in burton mesh picked up from burton_mesh_traits_t.
  static constexpr size_t num_domains = mesh_traits_t::num_domains;

  //============================================================================
  // Define basic types.
  //============================================================================

  //! Type for burton mesh vertices.
  using vertex_t = burton_vertex_t<dimension>;

  //! Type for burton mesh edges.
  using edge_t = burton_edge_t<dimension>;

  //! Type for burton mesh faces.
  using face_t = burton_face_t<dimension>;

  //! Type for burton mesh cells.
  using cell_t = burton_cell_t<dimension>;

  //! Type for burton mesh corners.
  using corner_t = burton_corner_t<dimension>;

  //! Type for burton mesh wedges.
  using wedge_t = burton_wedge_t<dimension>;

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


////////////////////////////////////////////////////////////////////////////////
//! \class burton_mesh_types_t burton_types.h
//! \brief A collection of type information needed to specialize the flecsi
//!   low-level mesh infrastructure for ALE methods.
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_mesh_types_t<3> {

  //============================================================================
  // Define local traits to satisfy mesh_topology requirements.
  //============================================================================

  //! the mesh traites
  using mesh_traits_t = burton_mesh_traits_t<3>;

  //! The dimension of the burton mesh picked up from burton_mesh_traits_t.
  static constexpr size_t dimension = mesh_traits_t::dimension;

  //! The number of domains in burton mesh picked up from burton_mesh_traits_t.
  static constexpr size_t num_domains = mesh_traits_t::num_domains;

  //============================================================================
  // Define basic types.
  //============================================================================

  //! Type for burton mesh vertices.
  using vertex_t = burton_vertex_t<dimension>;

  //! Type for burton mesh edges.
  using edge_t = burton_edge_t<dimension>;

  //! Type for burton mesh cells.
  using face_t = burton_face_t<dimension>;

  //! Type for burton mesh cells.
  using cell_t = burton_cell_t<dimension>;

  //! Type for burton mesh corners.
  using corner_t = burton_corner_t<dimension>;

  //! Type for burton mesh wedges.
  using wedge_t = burton_wedge_t<dimension>;

  //============================================================================
  // Specify mesh parameterizations.
  //============================================================================

  //! Definitions of burton mesh entities and their domain.
  // clang-format off
  using entity_types =
      std::tuple<
        std::pair<flecsi::domain_<0>, vertex_t>,
        std::pair<flecsi::domain_<0>, edge_t>,
        std::pair<flecsi::domain_<0>, face_t>,
        std::pair<flecsi::domain_<0>, cell_t>,
        std::pair<flecsi::domain_<1>, wedge_t>,
        std::pair<flecsi::domain_<1>, corner_t>
      >;

  //! Connectivities are adjacencies of entities within a single domain.
  using connectivities =
    std::tuple<
      std::tuple<flecsi::domain_<0>, vertex_t, edge_t>,
      std::tuple<flecsi::domain_<0>, vertex_t, face_t>,
      std::tuple<flecsi::domain_<0>, vertex_t, cell_t>,
      std::tuple<flecsi::domain_<0>, edge_t, vertex_t>,
      std::tuple<flecsi::domain_<0>, edge_t, face_t>,
      std::tuple<flecsi::domain_<0>, edge_t, cell_t>,
      std::tuple<flecsi::domain_<0>, face_t, vertex_t>,
      std::tuple<flecsi::domain_<0>, face_t, edge_t>,
      std::tuple<flecsi::domain_<0>, face_t, cell_t>,
      std::tuple<flecsi::domain_<0>, cell_t, vertex_t>,
      std::tuple<flecsi::domain_<0>, cell_t, face_t>,
      std::tuple<flecsi::domain_<0>, cell_t, edge_t>
      >;

  //! Bindings are adjacencies of entities across two domains.
  using bindings =
      std::tuple<
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>,   cell_t, corner_t>,
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>,   face_t, corner_t>,
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>,   edge_t, corner_t>,
        std::tuple<flecsi::domain_<0>, flecsi::domain_<1>, vertex_t, corner_t>,
        std::tuple<flecsi::domain_<1>, flecsi::domain_<0>, corner_t,   cell_t>,
        std::tuple<flecsi::domain_<1>, flecsi::domain_<0>, corner_t,   face_t>,
        std::tuple<flecsi::domain_<1>, flecsi::domain_<0>, corner_t,   edge_t>,
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
