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

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

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
  burton_wedge_t() = default;

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
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, bool is_right );

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
  
  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! \brief Is this wedge on the boundary.
  //! \return true if on boundary.
  auto is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }
  
  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  //============================================================================
  // Private Data
  //============================================================================

 private:

  //! centroid of the outer facetn
  point_t facet_centroid_ = 0;
  //! the normal of the outer facet
  vector_t facet_normal_ = 0;
  //! the area of the outer facet
  real_t facet_area_ = 0;
  //! the entity flags
  bitfield_t flags_;


};

////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 2d wedge
////////////////////////////////////////////////////////////////////////////////
template<>
template< typename MESH_TOPOLOGY >
void burton_wedge_t<2>::update( const MESH_TOPOLOGY * mesh, bool is_right )
{
  using math::abs;
  auto vs = mesh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = mesh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  const auto & e = es.front()->midpoint();
  const auto & v = vs.front()->coordinates();
  if ( is_right )
    facet_normal_ = { e[1] - v[1], v[0] - e[0] };
  else
    facet_normal_ = { v[1] - e[1], e[0] - v[0] };
  facet_area_ = abs(facet_normal_);
  facet_normal_ /= facet_area_;
  facet_centroid_ = 0.5 * ( e + v );
  set_boundary( es.front()->is_boundary() );
}

////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 3d wedge
////////////////////////////////////////////////////////////////////////////////
template<>
template< typename MESH_TOPOLOGY >
void burton_wedge_t<3>::update(const  MESH_TOPOLOGY* mesh, bool is_right)
{
  using math::abs;
  auto vs = mesh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = mesh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  auto fs = mesh->template entities<face_t::dimension, domain, face_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  assert( fs.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
  auto f = fs.front()->midpoint();
  if ( is_right )
    facet_normal_ = geom::shapes::triangle<num_dimensions>::normal( v, f, e );
  else 
    facet_normal_ = geom::shapes::triangle<num_dimensions>::normal( v, e, f );
  facet_area_ = abs(facet_normal_);
  facet_normal_ /= facet_area_;
  facet_centroid_ = geom::shapes::triangle<num_dimensions>::centroid( v, f, e );
  set_boundary( fs.front()->is_boundary() );
}


} // namespace
} // namespace
} // namespace
