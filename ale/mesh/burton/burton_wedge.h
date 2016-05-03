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
#include "ale/geom/shapes/triangle.h"
#include "ale/mesh/burton/burton_vertex.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {


//! forward decares
template< std::size_t N >
class burton_corner_t;


////////////////////////////////////////////////////////////////////////////////
//! \class burton_wedge_t burton_entity_types.h
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//!
//! \tparam N The domain of the wedge.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_wedge_t {};

////////////////////////////////////////////////////////////////////////////////
//! \class burton_wedge_t burton_entity_types.h
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//!
//! \tparam N The domain of the wedge.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_wedge_t<2>
  : public flecsi::mesh_entity_t<2, burton_mesh_traits_t<2>::num_domains>
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

  //! Physics vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_cell_t<num_dimensions>;

  //! the base corner type
  using corner_t = burton_corner_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_wedge_t(mesh_topology_base_t & mesh) {};

  //! dissallow copying
  burton_wedge_t( burton_wedge_t & ) = delete;
  burton_wedge_t & operator=( burton_wedge_t & ) = delete;

  //! dissallow moving
  burton_wedge_t( burton_wedge_t && ) = delete;
  burton_wedge_t & operator=( burton_wedge_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! Set the vertex that a wedge has.
  void set_vertex(vertex_t * vertex) { vertex_ = vertex; }

  //! Set the edge that a wedge has.
  void set_edge(edge_t * edge) { edge_ = edge; }

  //! Set the face that a wedge has.
  void set_face(face_t * face) { }

  //! Set the cell that a wedge is in.
  void set_cell(cell_t * cell) { cell_ = cell; }

  //! Set the corner that a wedge is in.
  void set_corner(corner_t * corner) { corner_ = corner; }

  //! Get the vertex that a wedge has.
  const vertex_t * vertex() const { return vertex_; }

  //! Get the edge that a wedge has.
  const edge_t * edge() const { return edge_; }

  //! Get the face that a wedge has.
  const edge_t * face() const { return edge_; }

  //! Get the cell that a wedge is in.
  const cell_t * cell() const { return cell_; }

  //! Get the corner that a wedge is in.
  const corner_t * corner() const { return corner_; }

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  vector_t facet_normal_left() const
  {
    auto e = edge()->midpoint();
    auto v = vertex()->coordinates();
    return { v[1] - e[1], e[0] - v[0] };
  }
  vector_t facet_normal_right() const
  {
    auto e = edge()->midpoint();
    auto v = vertex()->coordinates();
    return { e[1] - v[1], v[0] - e[0] };
  }

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) { }

  //============================================================================
  // Private Data
  //============================================================================

 private:

  vertex_t * vertex_;
  edge_t * edge_;
  cell_t * cell_;
  corner_t * corner_;

}; // struct burton_wedge_t

////////////////////////////////////////////////////////////////////////////////
//! \class burton_wedge_t burton_entity_types.h
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//!
//! \tparam N The domain of the wedge.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_wedge_t<3>
  : public flecsi::mesh_entity_t<2, burton_mesh_traits_t<3>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::mesh_topology_base_t;
 
  //! the mesh traits
  using mesh_traits_t = burton_mesh_traits_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = mesh_traits_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = mesh_traits_t::dimension;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Physics vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_cell_t<num_dimensions>;

  //! the base corner type
  using corner_t = burton_corner_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_wedge_t(mesh_topology_base_t & mesh) {};

  //! dissallow copying
  burton_wedge_t( burton_wedge_t & ) = delete;
  burton_wedge_t & operator=( burton_wedge_t & ) = delete;

  //! dissallow moving
  burton_wedge_t( burton_wedge_t && ) = delete;
  burton_wedge_t & operator=( burton_wedge_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! Set the vertex that a wedge has.
  void set_vertex(vertex_t * vertex) { vertex_ = vertex; }

  //! Set the edge that a wedge has.
  void set_edge(edge_t * edge) { edge_ = edge; }

  //! Set the face that a wedge has.
  void set_face(face_t * face) { face_ = face; }

  //! Set the cell that a wedge is in.
  void set_cell(cell_t * cell) { cell_ = cell; }

  //! Set the corner that a wedge is in.
  void set_corner(corner_t * corner) { corner_ = corner; }

  //! Get the vertex that a wedge has.
  const vertex_t * vertex() const { return vertex_; }

  //! Get the edge that a wedge has.
  const edge_t * edge() const { return edge_; }

  //! Get the face that a wedge has.
  const face_t * face() const { return face_; }

  //! Get the cell that a wedge is in.
  const cell_t * cell() const { return cell_; }

  //! Get the corner that a wedge is in.
  const corner_t * corner() const { return corner_; }

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  vector_t facet_normal_left() const
  {
    auto e = edge()->midpoint();
    auto v = vertex()->coordinates();
    auto f = face()->centroid();
    return geom::triangle<num_dimensions>::normal( v, f, e );
  }
  vector_t facet_normal_right() const
  {
    auto e = edge()->midpoint();
    auto v = vertex()->coordinates();
    auto f = face()->centroid();
    return geom::triangle<num_dimensions>::normal( v, e, f );
  }

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) { }

  //============================================================================
  // Private Data
  //============================================================================

 private:

  vertex_t * vertex_;
  edge_t * edge_;
  face_t * face_;
  cell_t * cell_;
  corner_t * corner_;

}; // struct burton_wedge_t

} // namespace
} // namespace
