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
#include "flecsi/mesh/mesh_types.h"

#include "ale/geom/geometric_shapes.h"
#include "ale/mesh/burton/burton_mesh_traits.h"
#include "ale/utils/errors.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
// Expose some types
////////////////////////////////////////////////////////////////////////////////

//! some flexi flags
using flecsi::persistent;
using flecsi::flecsi_internal;
using flecsi::flecsi_user_space;

//! Type for storing instance of template specialized low level mesh.
template < std::size_t N >
using burton_mesh_topology_t = 
  flecsi::mesh_topology_t< burton_mesh_types_t<N> >;



////////////////////////////////////////////////////////////////////////////////
//! \brief The base for all entity types
////////////////////////////////////////////////////////////////////////////////
template< std::size_t NUM_DIMS, std::size_t DIM >
class burton_entity_base_t : 
    public flecsi::mesh_entity_t<DIM, burton_mesh_traits_t<NUM_DIMS>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::mesh_topology_base_t;
 
  //! the flecsi mesh topology type
  using mesh_topology_t = burton_mesh_topology_t<NUM_DIMS>;

  //! the mesh traits
  using mesh_traits_t = burton_mesh_traits_t<NUM_DIMS>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = mesh_traits_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = mesh_traits_t::dimension;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_entity_base_t(mesh_topology_base_t & mesh) : mesh_(&mesh) {}

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //! get the mesh
  auto mesh() const
  { 
    return static_cast<const mesh_topology_t *>(mesh_); 
  }

  //============================================================================
  // Private Data
  //============================================================================

 private:

  //! a reference to the mesh topology
  const mesh_topology_base_t * mesh_ = nullptr;
};

////////////////////////////////////////////////////////////////////////////////
//! \class burton_vertex_t burton_entity_types.h
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//!
//! \tparam N The domain of the vertex.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_vertex_t : 
    public burton_entity_base_t<N, 0>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the entity base type
  using entity_base_t = burton_entity_base_t<N, 0>;

  //! the burton mesh topology type
  using typename entity_base_t::mesh_topology_base_t;

  //! the mesh traits
  using typename entity_base_t::mesh_traits_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! the bitfield type
  using bitfield_t = typename mesh_traits_t::bitfield_t;

  //! Number of domains in the burton mesh.
  using entity_base_t::num_domains;

  //! Number of domains in the burton mesh.
  using entity_base_t::num_dimensions;


  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_vertex_t(mesh_topology_base_t & mesh) : entity_base_t(mesh) 
  {}

  //! dissallow copying
  burton_vertex_t( burton_vertex_t & ) = delete;
  burton_vertex_t & operator=( burton_vertex_t & ) = delete;

  //! dissallow moving
  burton_vertex_t( burton_vertex_t && ) = delete;
  burton_vertex_t & operator=( burton_vertex_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  // using the mesh accessor from the base class
  using entity_base_t::mesh;

  //! \brief Set the coordinates for a vertex.
  //! \param[in] coordinates Coordinates value to set for vertex.
  void set_coordinates(const point_t & coordinates)
  {
    using flecsi::mesh_entity_base_t;
    auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
      "coordinates", mesh()->runtime_id() );
    c[mesh_entity_base_t<num_domains>::template id<0>()] = coordinates;
  } // set_coordinates

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  const point_t & coordinates() const
  {
    using flecsi::mesh_entity_base_t;
    const auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
      "coordinates", mesh()->runtime_id() );
    return c[mesh_entity_base_t<num_domains>::template id<0>()];
  } // coordinates

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  //! \remark this is the non const version
  point_t & coordinates()
  {
    using flecsi::mesh_entity_base_t;
    auto c = data_t::instance().template dense_accessor<point_t, flecsi_internal>(
      "coordinates", mesh()->runtime_id() );
    return c[mesh_entity_base_t<num_domains>::template id<0>()];
  } // coordinates

  //! is this a boundary
  bool is_boundary() const
  {
    using flecsi::mesh_entity_base_t;
    auto flag =
      data_t::instance().template dense_accessor<bitfield_t, flecsi_internal>(
        "point_flags", mesh()->runtime_id() );
    return flag[mesh_entity_base_t<num_domains>::template id<0>()].anybitset();
  }

}; // class burton_vertex_t

////////////////////////////////////////////////////////////////////////////////
//! \class burton_edge_t burton_entity_types.h
//! \brief The burton_edge_t type provides an interface for managing
//!   geometry and state associated with mesh edges.
//!
//! \tparam N The domain of the edge.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_edge_t :
    public burton_entity_base_t<N, 1>
{

  //============================================================================
  // Typedefs
  //============================================================================

  //! the entity base type
  using entity_base_t = burton_entity_base_t<N, 1>;

  //! the burton mesh topology type
  using typename entity_base_t::mesh_topology_base_t;

  //! the mesh traits
  using typename entity_base_t::mesh_traits_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! the bitfield type
  using bitfield_t = typename mesh_traits_t::bitfield_t;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;


  //! Number of domains in the burton mesh.
  using entity_base_t::num_domains;

  //! Number of domains in the burton mesh.
  using entity_base_t::num_dimensions;

  //============================================================================
  // Constructors
  //============================================================================

  //! the constructor
  burton_edge_t(mesh_topology_base_t & mesh) : entity_base_t(mesh) {}

  //! dissallow copying
  burton_edge_t( burton_edge_t & ) = delete;
  burton_edge_t & operator=( burton_edge_t & ) = delete;

  //! dissallow moving
  burton_edge_t( burton_edge_t && ) = delete;
  burton_edge_t & operator=( burton_edge_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! get the mesh
  using entity_base_t::mesh;

  //! the list of vertices
  auto vertices() const
  {
    auto vs = mesh()->template entities<0,0>(this);
    return vs;
  }

  //! the list of actual coordinates
  auto coordinates() const
  {
    auto vs = vertices();
    std::vector< point_t > coords;
    coords.reserve( vs.size() );
    for ( auto v : vs ) coords.emplace_back( v->coordinates() );
    return coords;
  }

  
  //! the edge midpoint
  point_t midpoint() const
  {
    auto vs = mesh()->template entities<0,0>(this);    
    return {0.5*(vs[0]->coordinates() + vs[1]->coordinates())};
  }

  //! the edge length
  real_t  length() const
  {
    using math::sqr;
    auto vs = mesh()->template entities<0,0>(this);
    auto & a = vs[0]->coordinates();
    auto & b = vs[1]->coordinates();
    return std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
  }

  //! in 2d, this doubles as a face, so the area is the same as the length
  //! \remark this is only enabled in 2d
  real_t area() const
  {
    return length();
  }

  //! the edge normal
  vector_t normal() const
  {
    auto vs = mesh()->template entities<0,0>(this);
    auto & a = vs[0]->coordinates();
    auto & b = vs[1]->coordinates();
    return { a[1] - b[1], b[0] - a[0] };
  }


  //! is this a boundary
  bool is_boundary() const
  {
    using flecsi::mesh_entity_base_t;
    auto flag =
      data_t::instance().template dense_accessor<bitfield_t, flecsi_internal>(
        "edge_flags", mesh()->runtime_id() );
    return flag[mesh_entity_base_t<num_domains>::template id<0>()].anybitset();
  }

}; // struct burton_edge_t

////////////////////////////////////////////////////////////////////////////////
//! \class burton_element_t burton_entity_types.h
//! \brief The burton_element_t type provides an interelement for managing and
//!   geometry and state associated with mesh elements.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t NUM_DIMS, std::size_t DIM  >
struct burton_element_t
  : public burton_entity_base_t<NUM_DIMS, DIM>
{

  //============================================================================
  // Typedefs
  //============================================================================

  //! the entity base type
  using entity_base_t = burton_entity_base_t<NUM_DIMS, DIM>;

  //! the burton mesh topology type
  using typename entity_base_t::mesh_topology_base_t;

  //! the mesh traits
  using typename entity_base_t::mesh_traits_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  // the flecsi id type
  using id_t = flecsi::id_t;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Number of domains in the burton mesh.
  using entity_base_t::num_domains;

  //! Number of domains in the burton mesh.
  using entity_base_t::num_dimensions;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(mesh_topology_base_t & mesh) : entity_base_t(mesh) 
  {};

  //! Destructor
  virtual ~burton_element_t() {}

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! get the mesh
  using entity_base_t::mesh;

  //! the list of vertices
  auto vertices() const
  {
    auto vs = mesh()->template entities<0,0>(this);
    return vs;
  }

  //! the list of actual coordinates
  auto coordinates() const
  {
    auto vs = vertices();
    std::vector< point_t > coords;
    coords.reserve( vs.size() );
    for ( auto v : vs ) coords.emplace_back( v->coordinates() );
    return coords;
  }


  //! the list of edges
  auto edges() const
  {
    auto es = mesh()->template entities<1,0>(this);
    return es;
  }

  //! the list of faces
  auto faces() const
  {
    auto fs = mesh()->template entities<DIM-1,0>(this);
    return fs;
  }

  //! the centroid
  virtual point_t centroid() const 
  { raise_runtime_error("you should never get here"); };


  //! the area of the element
  virtual real_t area() const
  { raise_runtime_error("you should never get here"); };

  //! the area of the element
  virtual real_t volume() const
  { area(); };

  //! the minimum length in the element
  virtual real_t min_length() const
  {
    // get the vertices
    auto vs = vertices();
    // get one of the edges as a reference
    auto es = edges();
    auto min_length = es.front()->length();
    // check each vertex combination
    for ( auto vi : vs ) {
      auto pi = vi->coordinates();
      for ( auto vj : vs ) {
        if ( vi == vj ) continue;
        auto pj = vj->coordinates();
        auto delta = pi - pj;
        min_length = std::min( abs(delta), min_length );
      }
    }
    // return the result
    return min_length;
  }


  //! the element type
  virtual geom::geometric_shapes_t type() const
  { raise_runtime_error("you should never get here"); };


  //! get the region id
  auto & region()
  {
    using flecsi::mesh_entity_base_t;
    auto cell_regions =
      data_t::instance().template dense_accessor<size_t, flecsi_internal>(
        "cell_region", mesh()->runtime_id() );
    return cell_regions[mesh_entity_base_t<num_domains>::template id<0>()];
  }

  //! get the region id
  auto region() const
  {
    using flecsi::mesh_entity_base_t;
    auto cell_regions =
      data_t::instance().template dense_accessor<size_t, flecsi_internal>(
        "cell_region", mesh()->runtime_id() );
    return cell_regions[mesh_entity_base_t<num_domains>::template id<0>()];
  }

  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[out] e Vector to fill with ids of the vertices making the entity.
  //! \param[in] v Vertex ids for the element.
  //! \param[in] vertex_count The number of vertices making up the entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  virtual std::vector<id_t> create_entities(
    size_t dim, id_t * e, id_t * v, size_t vertex_count) 
  { raise_runtime_error("you should never get here"); };

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] dim The topological dimension of the entity being bound.
  //! \param[in] ent_ids The entity ids of the entities making up the binding.
  //! \param[out] c The collection of the ids making up the bound entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  virtual std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, id_t ** ent_ids, 
    size_t * ent_counts, id_t * c ) 
  { raise_runtime_error("you should never get here"); };

}; // class burton_element_t


////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_face_t type provides an interface for managing and
//!   geometry and state associated with mesh faces.
//!
//! \tparam N The domain of the face.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_face_t = burton_element_t<N,N-1>;

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_cell_t type provides an interface for managing and
//!   geometry and state associated with mesh cells.
//!
//! \tparam N The domain of the cell.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_cell_t = burton_element_t<N,N>;


////////////////////////////////////////////////////////////////////////////////
//! \class burton_wedge_t burton_entity_types.h
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//!
//! \tparam N The domain of the wedge.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_wedge_t
  : public burton_entity_base_t<N, 2>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the entity base type
  using entity_base_t = burton_entity_base_t<N, 2>;

  //! the burton mesh topology type
  using typename entity_base_t::mesh_topology_base_t;

  //! the mesh traits
  using typename entity_base_t::mesh_traits_t;

  //! Physics vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<N>;

  //! the base edge type
  using edge_t = burton_edge_t<N>;

  //! the base cell type
  using cell_t = burton_cell_t<N>;

  //! the base corner type
  using corner_t = burton_corner_t<N>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_wedge_t(mesh_topology_base_t & mesh) : entity_base_t(mesh)
  {}

  //! dissallow copying
  burton_wedge_t( burton_wedge_t & ) = delete;
  burton_wedge_t & operator=( burton_wedge_t & ) = delete;

  //! dissallow moving
  burton_wedge_t( burton_wedge_t && ) = delete;
  burton_wedge_t & operator=( burton_wedge_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! Set the cell that a wedge is in.
  void set_cell(cell_t * cell) { cell_ = cell; }

  //! Set the edge that a wedge has.
  void set_edge(edge_t * edge) { edge_ = edge; }

  //! Set the vertex that a wedge has.
  void set_vertex(vertex_t * vertex) { vertex_ = vertex; }

  //! Set the corner that a wedge is in.
  void set_corner(corner_t * corner) { corner_ = corner; }

  //! Get the cell that a wedge is in.
  const cell_t * cell() const { return cell_; }

  //! Get the edge that a wedge has.
  const edge_t * edge() const { return edge_; }

  //! Get the vertex that a wedge has.
  const vertex_t * vertex() const { return vertex_; }

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

  //============================================================================
  // Private Data
  //============================================================================

 private:

  cell_t * cell_;
  edge_t * edge_;
  vertex_t * vertex_;
  corner_t * corner_;

}; // struct burton_wedge_t

////////////////////////////////////////////////////////////////////////////////
//! \class burton_corner_t burton_entity_types.h
//! \brief The burton_corner_t type provides an interface for managing and
//!   geometry and state associated with mesh corners.
//!
//! \tparam N The domain of the corner.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_corner_t
  : public burton_entity_base_t<N, 1>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the entity base type
  using entity_base_t = burton_entity_base_t<N, 1>;

  //! the burton mesh topology type
  using typename entity_base_t::mesh_topology_base_t;

  //! the mesh traits
  using typename entity_base_t::mesh_traits_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type containing coordinates of a vertex.
  using point_t = typename mesh_traits_t::point_t;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<N>;

  //! the base edge type
  using edge_t = burton_edge_t<N>;

  //! the base cell type
  using cell_t = burton_cell_t<N>;

  //! the base wedge type
  using wedge_t = burton_wedge_t<N>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_corner_t(mesh_topology_base_t & mesh) : entity_base_t(mesh) 
  {}

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

  //! Set the an edge that a corner has.
  void add_edge(edge_t * edge) { edges_.add(edge); }

  //! Set the vertex that a corner has.
  void set_vertex(vertex_t * vertex) { vertex_ = vertex; }

  //! Get the cell that a corner is in.
  const cell_t * cell() const { return cell_; }

  //! Get edge1 that a corner has.
  auto & edges() const { return edges_; }

  //! Get the vertex that a corner has.
  const vertex_t * vertex() const { return vertex_; }

  //! \brief Get the wedges for the mesh.
  //! \return The wedges in the mesh.
  auto & wedges() { return wedges_; } // wedges

  //============================================================================
  // Private Data
  //============================================================================

private:


  flecsi::entity_group<wedge_t> wedges_;
  flecsi::entity_group<edge_t> edges_;
  cell_t * cell_;
  vertex_t * vertex_;

}; // class burton_corner_t

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
