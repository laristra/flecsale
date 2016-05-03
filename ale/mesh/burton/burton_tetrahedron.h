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
#include "ale/geom/shapes/tetrahedron.h"
#include "ale/mesh/burton/burton_element.h"


namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \class burton_tetrahedron_t burton_entity_types.h
//!
//! \brief The burton_tetrahedron_t type provides a derived instance of
//!   burton_cell_t for 2D tetrahedron cells.
////////////////////////////////////////////////////////////////////////////////
class burton_tetrahedron_t : public burton_element_t<3,3>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the base cell type
  using base_t = burton_element_t<3,3>;

  //============================================================================
  // Constructors
  //============================================================================

  //! main constructor
  burton_tetrahedron_t(mesh_topology_base_t & mesh) : base_t(mesh)
  { }

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  point_t centroid() const override;

  //! the area of the cell
  real_t volume() const override;

  //! the cell type
  geom::geometric_shapes_t type() const override 
  { return geom::tetrahedron::shape; };

  //----------------------------------------------------------------------------
  //! \brief create_entities function for burton_tetrahedron_cell_t.
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_entities(
      size_t dim, id_t * e, id_t * v, size_t vertex_count)  override
  {

    assert( vertex_count == 4 );

    switch (dim) {

      //------------------------------------------------------------------------
      // Edges
    case (1): 
      // bottom
      e[0] = v[0];
      e[1] = v[1];

      e[2] = v[1];
      e[3] = v[2];

      e[4] = v[2];
      e[5] = v[0];

      // top
      e[6] = v[0];
      e[7] = v[3];

      e[8] = v[1];
      e[9] = v[3];

      e[10] = v[2];
      e[11] = v[3];

      return {2, 2, 2, 2, 2, 2};

      //------------------------------------------------------------------------
      // Faces
    case (2): 

      e[0]  = v[1];
      e[1]  = v[3];
      e[2]  = v[2];

      e[3]  = v[2];
      e[4]  = v[3];
      e[5]  = v[0];

      e[6]  = v[0];
      e[7]  = v[3];
      e[8]  = v[1];

      e[9]  = v[0];
      e[10] = v[1];
      e[11] = v[2];

      return {3, 3, 3, 3};
      
      //------------------------------------------------------------------------
      // Failure      
    default: 
      raise_runtime_error("Unknown entity type");

    } // switch


  } // create_entities

  //----------------------------------------------------------------------------
  /*!
    \brief create_bound_entities function for burton_tetrahedron_cell_t.

    \verbatim

    The following shows the labeling of the primitives making up a cell. Given
    vertices v*, edges e*, and center vertex cv.

    v3------e2-------v2
    |                 |
    |                 |
    |                 |
    |                 |
    e3      cv       e1
    |                 |
    |                 |
    |                 |
    |                 |
    v0------e0-------v1

    A wedge is defined by a vertex, an edge, and the cell itself. The wedge
    indexing is shown below.

    v3------e2-------v2
    | \      |      / |
    |   \  w6|w5  /   |
    |  w7 \  |  / w4  |
    |       \|/       |
    e3------cv-------e1
    |       /|\       |
    |  w0 /  |  \ w3  |
    |   /  w1|w2  \   |
    | /      |      \ |
    v0------e0-------v1

    A corner is defined by a vertex and two edges.

    c0 = {v0, e0, e3}
    c1 = {v1, e0, e1}
    c2 = {v2, e1, e2}
    c3 = {v3, e2, e3}

    \endverbatim
   */
  //----------------------------------------------------------------------------
  inline std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, id_t ** ent_ids, 
    size_t * ent_counts, id_t * c ) override
  {
    assert( ent_counts[0] == 4 );

    switch (dim) {
      // Corners
      // The right edge is always first
      case 1:
        // corner 0
        c[0] = ent_ids[0][0]; // vertex 0
        c[1] = ent_ids[1][0]; // edge 0, abuts vertex 0
        c[2] = ent_ids[1][3]; // edge 3, abuts vertex 0

        // corner 1
        c[3] = ent_ids[0][1]; // vertex 1
        c[4] = ent_ids[1][1]; // edge 1, abuts vertex 1
        c[5] = ent_ids[1][0]; // edge 0, abuts vertex 1

        // corner 2
        c[6] = ent_ids[0][2]; // vertex 2
        c[7] = ent_ids[1][2]; // edge 2, abuts vertex 2
        c[8] = ent_ids[1][1]; // edge 1, abuts vertex 2

        // corner 3
        c[9] = ent_ids[0][3]; // vertex 3
        c[10] = ent_ids[1][3]; // edge 3, abuts vertex 3
        c[11] = ent_ids[1][2]; // edge 2, abuts vertex 3

        return {3, 3, 3, 3};

      // Wedges
      case 2:

        // wedge 0
        c[0] = ent_ids[0][0]; // vertex 0
        c[1] = ent_ids[1][3]; // edge 3

        // wedge 1
        c[2] = ent_ids[0][0]; // vertex 0
        c[3] = ent_ids[1][0]; // edge 0

        // wedge 2
        c[4] = ent_ids[0][1]; // vertex 1
        c[5] = ent_ids[1][0]; // edge 0

        // wedge 3
        c[6] = ent_ids[0][1]; // vertex 1
        c[7] = ent_ids[1][1]; // edge 1

        // wedge 4
        c[8] = ent_ids[0][2]; // vertex 2
        c[9] = ent_ids[1][1]; // edge 1

        // wedge 5
        c[10] = ent_ids[0][2]; // vertex 2
        c[11] = ent_ids[1][2]; // edge 2

        // wedge 6
        c[12] = ent_ids[0][3]; // vertex 3
        c[13] = ent_ids[1][2]; // edge 2

        // wedge 7
        c[14] = ent_ids[0][3]; // vertex 3
        c[15] = ent_ids[1][3]; // edge 3

        return {2, 2, 2, 2, 2, 2, 2, 2};

      default:
        raise_runtime_error("Unknown bound entity type");
    } // switch
  } // create_bound_entities



};



} // namespace
} // namespace
