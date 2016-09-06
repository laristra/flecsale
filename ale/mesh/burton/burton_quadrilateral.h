/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file 
/// \brief Provides various quadrilateral element types. 
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "ale/geom/shapes/quadrilateral.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \brief Provides a derived instance of burton_cell_t for multidimensional 
//!        quadrilateral elements.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_quadrilateral_t {};

////////////////////////////////////////////////////////////////////////////////
//! \brief Provides a two-dimensional quadrilateral cell.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_quadrilateral_t<2> : public burton_element_t<2,2>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<2,2>;

  //============================================================================
  // Constructors
  //============================================================================

  //! main constructor
  burton_quadrilateral_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  point_t centroid() const override;

  //! the midpoint
  point_t midpoint() const override;

  //! the area of the cell
  real_t area() const override;

  //! the minimum length in the cell
  real_t min_length() const override;

  //! the cell type
  geom::shapes::geometric_shapes_t type() const override 
  { return geom::shapes::quadrilateral<num_dimensions>::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_quadrilateral_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities ) override
  {

    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, vertex_t::dimension );

    assert( v.size() == 4 );

    entities[0] = v[0];
    entities[1] = v[1];

    entities[2] = v[1];
    entities[3] = v[2];

    entities[4] = v[2];
    entities[5] = v[3];

    entities[6] = v[3];
    entities[7] = v[0];

    return {2, 2, 2, 2};
  } // create_entities

  //----------------------------------------------------------------------------
  //
  //  \brief create_bound_entities function for burton_quadrilateral_cell_t.
  //  
  //  \verbatim
  //  
  //  The following shows the labeling of the primitives making up a cell. Given
  //  vertices v*, edges e*, and center vertex cv.
  //  
  //  v3------e2-------v2
  //  |                 |
  //  |                 |
  //  |                 |
  //  |                 |
  //  e3      cv       e1
  //  |                 |
  //  |                 |
  //  |                 |
  //  |                 |
  //  v0------e0-------v1
  //  
  //  A wedge is defined by a vertex, an edge, and the cell itself. The wedge
  //  indexing is shown below.
  //  
  //  v3------e2-------v2
  //  | \      |      / |
  //  |   \  w6|w5  /   |
  //  |  w7 \  |  / w4  |
  //  |       \|/       |
  //  e3------cv-------e1
  //  |       /|\       |
  //  |  w0 /  |  \ w3  |
  //  |   /  w1|w2  \   |
  //  | /      |      \ |
  //  v0------e0-------v1
  //  
  //  A corner is defined by a vertex and two edges.
  //  
  //  c0 = {v0, e0, e3}
  //  c1 = {v1, e0, e1}
  //  c2 = {v2, e1, e2}
  //  c3 = {v3, e2, e3}
  //  
  //  \endverbatim
  //  
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * c )  override
  {

    auto verts = primal_conn.get_entity_vec( cell, vertex_t::dimension );
    auto edges = primal_conn.get_entity_vec( cell,   edge_t::dimension );

    assert( verts.size() == 4 );

    size_t i = 0;
    switch (dim) {

      //------------------------------------------------------------------------
      // Corners
      // The right edge is always first
    case 0:

      // corner 0
      c[i++] = verts[0]; // vertex 0
      c[i++] = edges[0]; // edge 0, abuts vertex 0
      c[i++] = edges[3]; // edge 3, abuts vertex 0
      // corner 1
      c[i++] = verts[1]; // vertex 1
      c[i++] = edges[1]; // edge 1, abuts vertex 1
      c[i++] = edges[0]; // edge 0, abuts vertex 1
      // corner 2
      c[i++] = verts[2]; // vertex 2
      c[i++] = edges[2]; // edge 2, abuts vertex 2
      c[i++] = edges[1]; // edge 1, abuts vertex 2
      // corner 3
      c[i++] = verts[3]; // vertex 3
      c[i++] = edges[3]; // edge 3, abuts vertex 3
      c[i++] = edges[2]; // edge 2, abuts vertex 3

      return {3, 3, 3, 3};

      //------------------------------------------------------------------------
      // Wedges
      // The right wedge is always first
    case 1: {

      auto corners = domain_conn.get_entity_vec( cell, corner_t::dimension );

      // wedge 0
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[0]; // edge 0, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 1
      c[i++] =   verts[0]; // vertex 0
      c[i++] =   edges[3]; // edge 3, abuts vertex 0
      c[i++] = corners[0]; // corner 0
      // wedge 2
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[1]; // edge 1, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 3
      c[i++] =   verts[1]; // vertex 1
      c[i++] =   edges[0]; // edge 0, abuts vertex 1
      c[i++] = corners[1]; // corner 1
      // wedge 4
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[2]; // edge 2, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 5
      c[i++] =   verts[2]; // vertex 2
      c[i++] =   edges[1]; // edge 1, abuts vertex 2
      c[i++] = corners[2]; // corner 2
      // wedge 6
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[3]; // edge 3, abuts vertex 3
      c[i++] = corners[3]; // corner 3
      // wedge 7
      c[i++] =   verts[3]; // vertex 3
      c[i++] =   edges[2]; // edge 2, abuts vertex 3
      c[i++] = corners[3]; // corner 3

      return {3, 3, 3, 3, 3, 3, 3, 3};
    }
      //------------------------------------------------------------------------
      // Failure
    default:
      raise_runtime_error("Unknown bound entity type");
    } // switch
  } // create_bound_entities


};



////////////////////////////////////////////////////////////////////////////////
//! \brief Provides a three dimensional quadrilateral face.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_quadrilateral_t<3> : public burton_element_t<3,2>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<3,2>;

  //============================================================================
  // Constructors
  //============================================================================

  //! main constructor
  burton_quadrilateral_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  point_t centroid() const override;

  //! the midpoint
  point_t midpoint() const override;

  //! the normal
  vector_t normal() const override;

  //! the area of the cell
  real_t area() const override;

  //! the minimum length in the cell
  real_t min_length() const override;

  //! the cell type
  geom::shapes::geometric_shapes_t type() const override 
  { return geom::shapes::quadrilateral<num_dimensions>::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_quadrilateral_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities ) override
  {

    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, vertex_t::dimension );

    assert( v.size() == 4 );

    entities[0] = v[0];
    entities[1] = v[1];

    entities[2] = v[1];
    entities[3] = v[2];

    entities[4] = v[2];
    entities[5] = v[3];

    entities[6] = v[3];
    entities[7] = v[0];

    return {2, 2, 2, 2};

  } // create_entities


};



} // namespace
} // namespace
