/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Provides the wedge type for use in burton_mesh_t.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "ale/geom/shapes/triangle.h"
#include "ale/mesh/burton/burton_vertex.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {


//! forward decares
template< std::size_t N >
class burton_corner_t;


////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//!
//! \tparam N The dimension of the wedge.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_wedge_t {};

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//! \remark This is a two-dimensional specialization.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_wedge_t<2>
  : public flecsi::topology::mesh_entity_t<1, burton_config_t<2>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::topology::mesh_topology_base_t;
 
  //! the mesh traits
  using config_t = burton_config_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 1;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //! Coordinate point type.
  using point_t = typename config_t::point_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_cell_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  // default constructor
  burton_wedge_t(mesh_topology_base_t & mesh) : mesh_(&mesh) {};

  // dissallow copying
  burton_wedge_t( burton_wedge_t & ) = delete;
  burton_wedge_t & operator=( burton_wedge_t & ) = delete;

  // dissallow moving
  burton_wedge_t( burton_wedge_t && ) = delete;
  burton_wedge_t & operator=( burton_wedge_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================


  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  vector_t facet_normal_left() const;
  //! \copydoc facet_normal_left
  vector_t facet_normal_right() const;
  
  //! \brief The "left" orientation of the facet normal.
  //! \return The facet normal vector.
  static vector_t facet_normal_left(
    const point_t & v, //!< [in] The associated vertex coordinate
    const point_t & e, //!< [in] The associated edge midpoint
    const point_t &    //!< [in] The associated face midpoint (not used in 2d)
  ) {
    return { v[1] - e[1], e[0] - v[0] };
  }

  //! \brief The "right" orientation of the facet normal.
  //! \return The facet normal vector.
  static vector_t facet_normal_right(
    const point_t & v, //!< [in] The associated vertex coordinate
    const point_t & e, //!< [in] The associated edge midpoint
    const point_t &    //!< [in] The associated face midpoint (not used in 2d)
  ) {
    return { e[1] - v[1], v[0] - e[0] };
  }
  

  //! \brief Get the cell facet centroid
  //! \return Cell facet centroid.
  point_t facet_centroid() const;
  //! \brief Get the cell facet midpoint
  //! \return Cell facet midpoint.
  point_t facet_midpoint() const;

  //! \brief reset the mesh pointer
  //! \param [in] mesh The new mesh to point to.
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh;
  }


  //! \brief Is this wedge on the boundary.
  //! \return true if on boundary.
  bool is_boundary() const;

  //============================================================================
  // Private Data
  //============================================================================

 private:

  //! a reference to the mesh topology
  mesh_topology_base_t * mesh_ = nullptr;

};

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_wedge_t type provides an interface for managing and
//!   geometry and state associated with mesh wedges.
//! \remark This is a three-dimensional specialization.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_wedge_t<3>
  : public flecsi::topology::mesh_entity_t<1, burton_config_t<3>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the flecsi mesh topology type
  using mesh_topology_base_t =  flecsi::topology::mesh_topology_base_t;
 
  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 1;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //! Coordinate point type.
  using point_t = typename config_t::point_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_cell_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  // default constructor
  burton_wedge_t(mesh_topology_base_t & mesh) : mesh_(&mesh) {};

  // dissallow copying
  burton_wedge_t( burton_wedge_t & ) = delete;
  burton_wedge_t & operator=( burton_wedge_t & ) = delete;

  // dissallow moving
  burton_wedge_t( burton_wedge_t && ) = delete;
  burton_wedge_t & operator=( burton_wedge_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  vector_t facet_normal_left() const;
  //! \copydoc facet_normal_left
  vector_t facet_normal_right() const;

  //! \brief The "left" orientation of the facet normal.
  //! \return The facet normal vector.
  static vector_t facet_normal_left(
    const point_t & v, //!< [in] The associated vertex coordinate
    const point_t & e, //!< [in] The associated edge midpoint
    const point_t & f  //!< [in] The associated face midpoint
  ) {
    return geom::shapes::triangle<num_dimensions>::normal( v, e, f );
  }

  //! \brief The "right" orientation of the facet normal.
  //! \return The facet normal vector.
  static vector_t facet_normal_right(
    const point_t & v, //!< [in] The associated vertex coordinate
    const point_t & e, //!< [in] The associated edge midpoint
    const point_t & f  //!< [in] The associated face midpoint
  ) {
    return geom::shapes::triangle<num_dimensions>::normal( v, f, e );
  }

  //! \brief Get the cell facet centroid
  //! \return Cell facet centroid.
  point_t facet_centroid() const;

  //! \brief Get the cell facet midpoint
  //! \return Cell facet midpoint.
  point_t facet_midpoint() const;

  //! \brief reset the mesh pointer
  //! \param [in] mesh The new mesh to point to.
  void reset(mesh_topology_base_t & mesh) 
  { 
    mesh_ = &mesh;
  }

  //! \brief Is this wedge on the boundary.
  //! \return true if on boundary.
  bool is_boundary() const;

  //============================================================================
  // Private Data
  //============================================================================

 private:

  //! a reference to the mesh topology
  mesh_topology_base_t * mesh_ = nullptr;

}; // struct burton_wedge_t

} // namespace
} // namespace
