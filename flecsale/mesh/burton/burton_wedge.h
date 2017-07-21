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
#include "flecsale/geom/shapes/triangle.h"
#include "flecsale/mesh/burton/burton_vertex.h"
#include "flecsale/mesh/burton/burton_element.h"


namespace flecsale {
namespace mesh {
namespace burton {


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
class burton_wedge_t
  : public flecsi::topology::mesh_entity_t<1, burton_config_t<N>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<N>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 1;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage_t<num_dimensions, num_domains>;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base_t< mesh_storage_t >;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

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

  //! \brief update the wedge geometry
  //! \param [in] is_right  This wedge is the right orientation when true.
  void update( bool is_right );

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  const auto & facet_normal() const
  { return facet_normal_; }

  //! \brief the facet area
  auto facet_area() const
  { return facet_area_; }

  //! \brief Get the cell facet centroid
  const auto & facet_centroid() const
  { return facet_centroid_; }
  
  //! \brief Get the cell facet midpoint
  const auto & facet_midpoint() const
  { return facet_centroid_; }
  

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

  //! geometric information
  point_t facet_centroid_ = 0;
  vector_t facet_normal_ = 0;
  real_t facet_area_ = 0;

};

} // namespace
} // namespace
} // namespace
