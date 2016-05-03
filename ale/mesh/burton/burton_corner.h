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
 * \file burton_entity_types.h
 * \authors bergen
 * \date Initial file creation: Nov 15, 2015
 ******************************************************************************/

#pragma once

//! user includes
#include "ale/mesh/burton/burton_vertex.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {



//! forward decares
template< std::size_t N >
class burton_wedge_t;

////////////////////////////////////////////////////////////////////////////////
//! \class burton_corner_t burton_entity_types.h
//! \brief The burton_corner_t type provides an interface for managing and
//!   geometry and state associated with mesh corners.
//!
//! \tparam N The domain of the corner.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_corner_t
  : public flecsi::mesh_entity_t<1, burton_mesh_traits_t<N>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::mesh_topology_base_t;
 
  //! the mesh traits
  using mesh_traits_t = burton_mesh_traits_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = mesh_traits_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = mesh_traits_t::dimension;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type containing coordinates of a vertex.
  using point_t = typename mesh_traits_t::point_t;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<N>;

  //! the base cell type
  using cell_t = burton_cell_t<N>;

  //! the base wedge type
  using wedge_t = burton_wedge_t<N>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_corner_t(mesh_topology_base_t & mesh) {}

  //! dissallow copying
  burton_corner_t( burton_corner_t & ) = delete;
  burton_corner_t & operator=( burton_corner_t & ) = delete;

  //! dissallow moving
  burton_corner_t( burton_corner_t && ) = delete;
  burton_corner_t & operator=( burton_corner_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief Add a wedge to the mesh.
  //!
  //! \param[in] w The wedge to add to the mesh.
  void add_wedge(wedge_t * w)
  {
    wedges_.add(w);
    w->set_corner(this);
  }

  //! Set the cell that a corner is in.
  void set_cell(cell_t * cell) { cell_ = cell; }

  //! Set the vertex that a corner has.
  void set_vertex(vertex_t * vertex) { vertex_ = vertex; }

  //! Get the cell that a corner is in.
  const cell_t * cell() const { return cell_; }

  //! Get the vertex that a corner has.
  const vertex_t * vertex() const { return vertex_; }

  //! \brief Get the wedges for the mesh.
  //! \return The wedges in the mesh.
  auto & wedges() { return wedges_; } // wedges

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) { }

  //============================================================================
  // Private Data
  //============================================================================

private:


  flecsi::entity_group<wedge_t> wedges_;
  cell_t * cell_;
  vertex_t * vertex_;

}; // class burton_corner_t

} // namespace mesh
} // namespace ale
