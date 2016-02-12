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

#include "flecsi/mesh/mesh_types.h"
#include "ale/mesh/burton/burton_mesh_traits.h"

/*!
 * \file burton_entity_types.h
 * \authors bergen
 * \date Initial file creation: Nov 15, 2015
 */

namespace ale {
namespace mesh {

//! the flecsi mesh topology type
using flecsi::mesh_topology_t;
//! the flecsi mesh topology type
using mesh_topology_base_t = flecsi::mesh_topology_base_t;

//! the flexi domain type
using flecsi::domain_;
//! the flexi domain entity type
using flecsi::domain_entity;
//! the flecsi entity group
using flecsi::entity_group;


//! some flexi flags
using flecsi::persistent;
using flecsi::flecsi_internal;
using flecsi::flecsi_user_space;

// some flecsi types
using flecsi::id_t;
using flecsi::id_vector_t;

// forward declarations
class burton_corner_t;
class burton_wedge_t;

/*----------------------------------------------------------------------------*
 * class burton_vertex_t
 *----------------------------------------------------------------------------*/

/*!
  \class burton_vertex_t burton_entity_types.h
  \brief The burton_vertex_t type provides an interface for managing
    geometry and state associated with mesh vertices.

  \tparam N The domain of the vertex.
 */
class burton_vertex_t
  : public flecsi::mesh_entity_t<0, burton_mesh_traits_t::num_domains>
{
public:
  
  //! Type containing coordinates of the vertex.
  using point_t = burton_mesh_traits_t::point_t;

  //! Handle for accessing state at vertex.
  using state_t = burton_mesh_traits_t::mesh_state_t;

  //! Number of domains in the burton mesh.
  static constexpr size_t num_domains = burton_mesh_traits_t::num_domains;

  //! Constructor
  burton_vertex_t(state_t & state) : state_(state) {}

  /*!
    \brief Set the coordinates for a vertex.

    \param[in] coordinates Coordinates value to set for vertex.
   */
  void set_coordinates(const point_t & coordinates)
  {
    auto c = state_.dense_accessor<point_t, flecsi_internal>("coordinates");
    c[mesh_entity_base_t<num_domains>::template id<0>()] = coordinates;
  } // set_coordinates

  /*!
    \brief Get the coordinates at a vertex from the state handle.
    \return coordinates of vertex.
   */
  const point_t & coordinates() const
  {
    const auto c =
        state_.dense_accessor<point_t, flecsi_internal>("coordinates");
    return c[mesh_entity_base_t<num_domains>::template id<0>()];
  } // coordinates

private:

  state_t & state_;

}; // class burton_vertex_t

/*----------------------------------------------------------------------------*
 * class burton_edge_t
 *----------------------------------------------------------------------------*/

/*!
  \class burton_edge_t burton_entity_types.h
  \brief The burton_edge_t type provides an interface for managing
    geometry and state associated with mesh edges.

  \tparam N The domain of the edge.
 */
struct burton_edge_t
  : public flecsi::mesh_entity_t<1, burton_mesh_traits_t::num_domains> {
  
  //! Type containing coordinates of the vertex.
  using point_t = burton_mesh_traits_t::point_t;

  burton_edge_t(mesh_topology_base_t & mesh)
    : mesh_(mesh) {}

  point_t midpoint() const;

private:

  //! a reference to the mesh topology
  mesh_topology_base_t & mesh_;

}; // struct burton_edge_t

/*----------------------------------------------------------------------------*
 * class burton_cell_t
 *----------------------------------------------------------------------------*/

/*!
  \class burton_cell_t burton_entity_types.h
  \brief The burton_cell_t type provides an interface for managing and
    geometry and state associated with mesh cells.
 */
struct burton_cell_t
  : public flecsi::mesh_entity_t<2, burton_mesh_traits_t::num_domains> {


  //! Type containing coordinates of the vertex.
  using point_t = burton_mesh_traits_t::point_t;

  //! Constructor
  burton_cell_t() = default;


  //! Destructor
  virtual ~burton_cell_t() {}

  virtual point_t centroid() const = 0;

  /*!
    \brief create_entities is a function that creates entities
      of topological dimension dim, using vertices v, and puts the vertices
      in e. See, e.g., burton_quadrilateral_cell_t for an implementation of
      this pure virtual function.

    \param[in] dim The topological dimension of the entity to create.
    \param[out] e Vector to fill with ids of the vertices making the entity.
    \param[in] v Vertex ids for the cell.
    \param[in] vertex_count The number of vertices making up the entity.

    \return A pair with a) the number of vertex collections making up the
      entity and b) the number of vertices per collection.
   */
  virtual id_vector_t
    create_entities( size_t dim, 
                     id_t * e, 
                     id_t * v, 
                     size_t vertex_count) = 0;

  /*!
    \brief create_bound_entities binds mesh entities across domains.
      See, e.g., burton_quadrilateral_cell_t for an implementation of
      this pure virtual function.

    \param[in] dim The topological dimension of the entity being bound.
    \param[in] ent_ids The entity ids of the entities making up the binding.
    \param[out] c The collection of the ids making up the bound entity.

    \return A pair with a) the number of entity collections making up the
      binding and b) the number of entities per collection.
   */
  virtual id_vector_t
    create_bound_entities( size_t from_domain, 
                           size_t to_domain,
                           size_t dim, 
                           id_t ** ent_ids,
                           id_t * c) = 0;

}; // class burton_cell_t

/*----------------------------------------------------------------------------*
 * class burton_quadrilateral_t
 *----------------------------------------------------------------------------*/

/*!
  \class burton_quadrilateral_t burton_entity_types.h
  \brief The burton_quadrilateral_t type provides a derived instance of
    burton_cell_t for 2D quadrilateral cells.
 */
class burton_quadrilateral_cell_t : public burton_cell_t
{
public:

  burton_quadrilateral_cell_t(mesh_topology_base_t & mesh)
    : mesh_(mesh) {}
 
  point_t centroid() const override;

  /*!
    \brief create_entities function for burton_quadrilateral_cell_t.
   */
  id_vector_t 
    create_entities( size_t dim, 
                     id_t * e, 
                     id_t * v, 
                     size_t vertex_count) override 
  {
    e[0] = v[0];
    e[1] = v[1];

    e[2] = v[1];
    e[3] = v[2];

    e[4] = v[2];
    e[5] = v[3];

    e[6] = v[3];
    e[7] = v[0];

    return {2, 2, 2, 2};
  } // create_entities
  
  /*!
    \brief create_bound_entities function for burton_quadrilateral_cell_t.
   */
  id_vector_t
  create_bound_entities( size_t from_domain, 
                         size_t to_domain, 
                         size_t dim,
                         id_t **ent_ids, 
                         id_t * c) override 
  {

    switch(dim) {
      // Corners
      case 1:
        // corner 0
        c[0] = ent_ids[0][0]; // vertex 0
        c[1] = ent_ids[1][0]; // edge 0
        c[2] = ent_ids[1][3]; // edge 3

        // corner 1
        c[4] = ent_ids[0][1]; // vertex 1
        c[5] = ent_ids[1][0]; // edge 0
        c[6] = ent_ids[1][1]; // edge 1

        // corner 2
        c[8] = ent_ids[0][2]; // vertex 2
        c[9] = ent_ids[1][1]; // edge 1
        c[10] = ent_ids[1][2]; // edge 2

        // corner 3
        c[12] = ent_ids[0][3]; // vertex 3
        c[13] = ent_ids[1][2]; // edge 2
        c[14] = ent_ids[1][3]; // edge 3

        return {4, 4, 4, 4};

#if 0 // Wedges are currently only referenced through corners
      // so this logic is unused for the time being...

      // Wedges
      case 2:
        c.resize(16);

        // wedge 0
        c[0] = ent_ids[0]; // vertex 0
        c[1] = ent_ids[7]; // edge 3

        // wedge 1
        c[2] = ent_ids[0]; // vertex 0
        c[3] = ent_ids[4]; // edge 0

        // wedge 2
        c[4] = ent_ids[1]; // vertex 1
        c[5] = ent_ids[4]; // edge 0

        // wedge 3
        c[6] = ent_ids[1]; // vertex 1
        c[7] = ent_ids[5]; // edge 1

        // wedge 4
        c[8] = ent_ids[2]; // vertex 2
        c[9] = ent_ids[5]; // edge 1

        // wedge 5
        c[10] = ent_ids[2]; // vertex 2
        c[11] = ent_ids[6]; // edge 2

        // wedge 6
        c[12] = ent_ids[3]; // vertex 3
        c[13] = ent_ids[6]; // edge 2

        // wedge 7
        c[14] = ent_ids[3]; // vertex 3
        c[15] = ent_ids[7]; // edge 3

        return {2, 2, 2, 2, 2, 2, 2, 2};
#endif

      default:
        assert(false && "Unknown bound entity type");
    } // switch
  } // create_bound_entities



private:

  //! a reference to the mesh topology
  mesh_topology_base_t & mesh_;

}; // class burton_quadrilateral_cell_t

/*----------------------------------------------------------------------------*
 * class burton_wedge_t
 *----------------------------------------------------------------------------*/

/*!
  \class burton_wedge_t burton_entity_types.h
  \brief The burton_wedge_t type provides an interface for managing and
    geometry and state associated with mesh wedges.

  \tparam N The domain of the wedge.
 */
class burton_wedge_t
  : public flecsi::mesh_entity_t<2, burton_mesh_traits_t::num_domains>
{
public:

  //! Physics vector type.
  using vector_t = burton_mesh_traits_t::vector_t;

  //! Set the corner that a wedge is in.
  void set_corner(burton_corner_t *corner) { corner_ = corner; }

  //! Get the corner that a wedge is in.
  burton_corner_t * corner() { return corner_; }

  /*!
    \brief Get the side facet normal for the wedge.
    \return Side facet normal vector.
   */
  vector_t side_facet_normal();

  /*!
    \brief Get the cell facet normal for the wedge.
    \return Cell facet normal vector.
   */
  vector_t cell_facet_normal();

private:

  burton_corner_t * corner_;

}; // struct burton_wedge_t

/*----------------------------------------------------------------------------*
 * class burton_corner_t
 *----------------------------------------------------------------------------*/

/*!
  \class burton_corner_t burton_entity_types.h
  \brief The burton_corner_t type provides an interface for managing and
    geometry and state associated with mesh corners.

  \tparam N The domain of the corner.
 */
class burton_corner_t
  : public flecsi::mesh_entity_t<1, burton_mesh_traits_t::num_domains>
{
public:

  burton_corner_t(mesh_topology_base_t & mesh)
    : mesh_(mesh) {}

  /*!
    \brief Add a wedge to the mesh.

    \param[in] w The wedge to add to the mesh.
   */
  void add_wedge(burton_wedge_t *w) {
    wedges_.add(w);
    w->set_corner(this);
  }

  /*!
    \brief Get the wedges for the mesh.
    \return The wedges in the mesh.
   */
  entity_group<burton_wedge_t> &wedges() { return wedges_; } // wedges

private:

  entity_group<burton_wedge_t> wedges_;

  //! a reference to the mesh topology
  mesh_topology_base_t & mesh_;

}; // class burton_corner_t

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
