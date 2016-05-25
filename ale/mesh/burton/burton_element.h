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

#pragma once


//! user includes
#include "ale/geom/shapes/geometric_shapes.h"
#include "ale/mesh/burton/burton_mesh_traits.h"
#include "ale/utils/errors.h"
#include "flecsi/mesh/mesh_types.h"

namespace ale {
namespace mesh {

//! forward decares
template< std::size_t N >
class burton_corner_t;

////////////////////////////////////////////////////////////////////////////////
//! \class burton_edge_t burton_entity_types.h
//! \brief The burton_edge_t type provides an interface for managing
//!   geometry and state associated with mesh edges.
//!
//! \tparam N The domain of the edge.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t NUM_DIMS, std::size_t DIM  >
struct burton_element_t { };

////////////////////////////////////////////////////////////////////////////////
//! \class burton_edge_t burton_entity_types.h
//! \brief The burton_edge_t type provides an interface for managing
//!   geometry and state associated with mesh edges.
//!
//! \tparam N The domain of the edge.
//!
//! \remark this is a 2d edge
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<2,1> :
    public flecsi::mesh_entity_t<1, burton_mesh_traits_t<2>::num_domains>
{

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

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! The bitfield.
  using bitfield_t = typename mesh_traits_t::bitfield_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! the constructor
  burton_element_t(mesh_topology_base_t & mesh) : mesh_(&mesh) {}

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the list of vertices
  auto vertices() const;

  //! the list of actual coordinates
  point_list_t coordinates() const;
  
  //! the edge midpoint
  point_t midpoint() const;

  //! \brief in 2d, this doubles as a face, so the centroid is the same as
  //!    the midpoint
  //! \remark this is only enabled in 2d
  point_t centroid() const;

  //! the edge length
  real_t  length() const;

  //! in 2d, this doubles as a face, so the area is the same as the length
  //! \remark this is only enabled in 2d
  real_t area() const;

  //! the edge normal
  vector_t normal() const;

  //! is this a boundary
  bool is_boundary() const;

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  const mesh_topology_base_t * mesh_ = nullptr;

}; // struct burton_edge_t


////////////////////////////////////////////////////////////////////////////////
//! \class burton_edge_t burton_entity_types.h
//! \brief The burton_edge_t type provides an interface for managing
//!   geometry and state associated with mesh edges.
//!
//! \tparam N The domain of the edge.
//!
//! \remark this is a 3d edge
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,1> :
    public flecsi::mesh_entity_t<1, burton_mesh_traits_t<3>::num_domains>
{

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

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! a bitfield type
  using bitfield_t = typename mesh_traits_t::bitfield_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! the constructor
  burton_element_t(mesh_topology_base_t & mesh) : mesh_(&mesh) {}

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the list of actual coordinates
  point_list_t coordinates() const;
  
  //! the edge midpoint
  point_t midpoint() const;

  //! the edge length
  real_t  length() const;

  //! is this a boundary
  bool is_boundary() const;

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  const mesh_topology_base_t * mesh_ = nullptr;

}; // struct burton_edge_t

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_face_t type provides an interface for managing and
//!   geometry and state associated with mesh faces.
//!
//! \tparam N The domain of the face.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_edge_t = burton_element_t<N,1>;

////////////////////////////////////////////////////////////////////////////////
//! \class burton_element_t burton_entity_types.h
//! \brief The burton_element_t type provides an interelement for managing and
//!   geometry and state associated with mesh elements.
//!
//! \remark this is a 2d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<2,2>
  : public flecsi::mesh_entity_t<2, burton_mesh_traits_t<2>::num_domains>
{
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

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type of floating point.
  using size_t = typename mesh_traits_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  // the flecsi id type
  using id_t = flecsi::id_t;
  using connectivity_t = flecsi::domain_connectivity<num_dimensions>;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base corner type
  using corner_t = burton_corner_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(mesh_topology_base_t & mesh) : mesh_(&mesh) 
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

  //! the list of actual coordinates
  point_list_t coordinates() const;

  //! the centroid
  virtual point_t centroid() const 
  { raise_runtime_error("you should never get here"); };

  //! the edge midpoint
  virtual point_t midpoint() const
  { raise_runtime_error("you should never get here"); };

  //! the area of the element
  virtual real_t area() const
  { raise_runtime_error("you should never get here"); };

  //! the area of the element
  virtual real_t volume() const
  { area(); };

  //! the minimum length in the element
  virtual real_t min_length() const;

  //! set the region id
  size_t & region();

  //! get the region id
  size_t region() const;


  //! the element type
  virtual geom::geometric_shapes_t type() const
  { raise_runtime_error("you should never get here"); };


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
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities ) 
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
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities ) 
  { raise_runtime_error("you should never get here"); };


  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //! \brief reset the mesh pointer
  const mesh_topology_base_t * mesh() const
  { 
    return mesh_; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  const mesh_topology_base_t * mesh_ = nullptr;

}; // class burton_element_t


////////////////////////////////////////////////////////////////////////////////
//! \class burton_element_t burton_entity_types.h
//! \brief The burton_element_t type provides an interelement for managing and
//!   geometry and state associated with mesh elements.
//!
//! \remark this is a 3d face
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,2>
  : public flecsi::mesh_entity_t<2, burton_mesh_traits_t<3>::num_domains>
{
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

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type of floating point.
  using size_t = typename mesh_traits_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  // the flecsi id type
  using id_t = flecsi::id_t;
  using connectivity_t = flecsi::domain_connectivity<num_dimensions>;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(mesh_topology_base_t & mesh) : mesh_(&mesh) 
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

  //! the list of actual coordinates
  point_list_t coordinates( bool reverse = false ) const;

  //! the centroid
  virtual point_t centroid() const 
  { raise_runtime_error("you should never get here"); };

  //! the midpoint used in tesselating the element
  virtual point_t midpoint() const 
  { raise_runtime_error("you should never get here"); };

  //! the direction the face is oriented in
  char direction() const
  { return direction_; }

  //! flip the direction the face is oriented in
  void flip() 
  { direction_ = - direction_; }

  //! flip the direction the face is oriented in
  bool flipped() const
  { return direction_ < 0; }

  //! the normal
  virtual vector_t normal() const 
  { raise_runtime_error("you should never get here"); };

  //! the area of the element
  virtual real_t area() const
  { raise_runtime_error("you should never get here"); };

  //! the minimum length in the element
  virtual real_t min_length() const;

  //! the element type
  virtual geom::geometric_shapes_t type() const
  { raise_runtime_error("you should never get here"); };

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
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities ) 
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
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities )
  { raise_runtime_error("you should never get here"); };

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //! \brief reset the mesh pointer
  const mesh_topology_base_t * mesh() const
  { 
    return mesh_; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  const mesh_topology_base_t * mesh_ = nullptr;

  //! a flag to indicate the direction of the face
  char direction_ = 1;

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
//! \class burton_element_t burton_entity_types.h
//! \brief The burton_element_t type provides an interelement for managing and
//!   geometry and state associated with mesh elements.
//!
//! \remark this is a 3d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,3>
  : public flecsi::mesh_entity_t<3, burton_mesh_traits_t<3>::num_domains>
{
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

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Handle for accessing state at vertex.
  using data_t = typename mesh_traits_t::data_t;

  //! Type of floating point.
  using real_t = typename mesh_traits_t::real_t;

  //! Type of floating point.
  using size_t = typename mesh_traits_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename mesh_traits_t::point_t;
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  // the flecsi id type
  using id_t = flecsi::id_t;
  using connectivity_t = flecsi::domain_connectivity<num_dimensions>;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_element_t<3,3>;

  //! the base corner type
  using corner_t = burton_corner_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(mesh_topology_base_t & mesh) : mesh_(&mesh)
  {}

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

  //! the list of actual coordinates
  point_list_t coordinates() const;

  //! the centroid
  virtual point_t centroid() const 
  { raise_runtime_error("you should never get here"); };

  //! the midpoint used in tesselating the element
  virtual point_t midpoint() const 
  { raise_runtime_error("you should never get here"); };

  //! the area of the element
  virtual real_t volume() const
  { raise_runtime_error("you should never get here"); };

  //! the minimum length in the element
  virtual real_t min_length() const;


  //! set the region id
  size_t & region();

  //! get the region id
  size_t region() const;

  //! the element type
  virtual geom::geometric_shapes_t type() const
  { raise_runtime_error("you should never get here"); };


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
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities )
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
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities )
  { raise_runtime_error("you should never get here"); };
 

  //! \brief reset the mesh pointer
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh; 
  }

  //! \brief reset the mesh pointer
  const mesh_topology_base_t * mesh() const
  { 
    return mesh_; 
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! a reference to the mesh topology
  const mesh_topology_base_t * mesh_ = nullptr;

}; // class burton_element_t



////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_cell_t type provides an interface for managing and
//!   geometry and state associated with mesh cells.
//!
//! \tparam N The domain of the cell.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_cell_t = burton_element_t<N,N>;

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
