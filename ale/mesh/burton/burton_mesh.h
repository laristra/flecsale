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
 * \file burton_mesh.h
 * \authors bergen
 * \date Initial file creation: Sep 02, 2015
 ******************************************************************************/

#pragma once

//! system includes
#include <set>
#include <string>

//! user includes
#include "flecsi/data/data.h"
#include "flecsi/execution/task.h"

#include "../../mesh/burton/burton_types.h"
#include "../../utils/errors.h"

namespace ale {
namespace mesh {


////////////////////////////////////////////////////////////////////////////////
/// \class burton_mesh_t burton_mesh.h
/// \brief A specialization of the flecsi low-level mesh topology, state and
///   execution models.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_mesh_t
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh types
  using mesh_types_t = burton_mesh_types_t<N>;

  //! the mesh traits
  using mesh_traits_t = typename mesh_types_t::mesh_traits_t;

  //! Type for storing instance of template specialized low level mesh.
  using mesh_t = flecsi::mesh_topology_t< mesh_types_t >;

  //! a compile string type
  using const_string_t = typename mesh_traits_t::const_string_t;
  
  //! a bitfield type
  using bitfield_t = typename mesh_traits_t::bitfield_t;

  //! Integer data type.
  using integer_t = typename mesh_traits_t::integer_t;

  //! Floating point data type.
  using real_t = typename mesh_traits_t::real_t;

  //! Floating point data type.
  using size_t = typename mesh_traits_t::size_t;

  //! Point data type.
  using point_t = typename mesh_traits_t::point_t;

  //! Physics vector type.
  using vector_t = typename mesh_traits_t::vector_t;

  //! Vertex type.
  using vertex_t = typename mesh_types_t::vertex_t;

  //! Edge type.
  using edge_t = typename mesh_types_t::edge_t;

  //! Cell type.
  using cell_t = typename mesh_types_t::cell_t;

  //! Wedge type.
  using wedge_t = typename mesh_types_t::wedge_t;

  //! Corner type.
  using corner_t = typename mesh_types_t::corner_t;

  //! \brief The locations of different bits that we set as flags
  enum bits : size_t 
  {
    boundary = 0 //!< the boundary bit
  };



  //! Type defining the execution policy.
#ifndef MESH_EXECUTION_POLICY
  using mesh_execution_t = flecsi::execution_t<>;
#else
  using mesh_execution_t = flecsi::execution_t<MESH_EXECUTION_POLICY>;
#endif

  //! Type defining the data attachment sites on the burton mesh.
  using attachment_site_t = typename mesh_traits_t::attachment_site_t;

  //! the data type
  using data_t = typename mesh_traits_t::data_t;

  //============================================================================
  // Dense Accessors
  //============================================================================

  //! \brief Accessor type.
  //! \tparam T The type of the underlying data to access.
  template <typename T>
  using dense_accessor_t = typename data_t::template dense_accessor_t<T>;

  //! \brief Register state for the named variable at the given attachment
  //!   site with attributes.
  //!
  //! \tparam T The type of the underlying data to access.
  //!
  //! \param[in] key A name for the state variable, e.g., "density".
  //! \param[in] site The data attachement site where the state variable should
  //!   be defined.  Valid sites are defined in flecsi::burton_mesh_traits_t.
  //! \param[in] attributes A bitfield specifying various attributes of the state.
  //!
  //! \return An accessor to the newly registered state.
  template <typename T>
  decltype(auto) register_state_(
    const const_string_t && key,
    attachment_site_t site, 
    typename bitfield_t::field_type_t attributes = 0x0 )
  {
    data_t & data_ = data_t::instance();

    switch (site) {
      case attachment_site_t::vertices:
        return data_.template register_state<T>(
            key, num_vertices(), mesh_.runtime_id(), 
            attachment_site_t::vertices, attributes 
        );
        break;
      case attachment_site_t::edges:
        return data_.template register_state<T>(
            key, num_edges(), mesh_.runtime_id(), 
            attachment_site_t::edges, attributes
        );
        break;
      case attachment_site_t::cells:
        return data_.template register_state<T>(
            key, num_cells(), mesh_.runtime_id(), 
            attachment_site_t::cells, attributes 
        );
        break;
      case attachment_site_t::corners:
        return data_.template register_state<T>(
            key, num_corners(), mesh_.runtime_id(), 
            attachment_site_t::corners, attributes
        );
        break;
      case attachment_site_t::wedges:
        return data_.template register_state<T>(
            key, num_wedges(), mesh_.runtime_id(), 
            attachment_site_t::wedges, attributes
        );
        break;
      default:
        assert(false && "Error: invalid state registration site.");
    } // switch

  } // register_state_

  //! \brief Access state associated with \e key.
  //!
  //! \tparam T Data type of underlying state.
  //! \tparam NS The namespace in which the state variable is registered.
  //!   See \ref data_t::register_state for additional information.
  //!
  //! \param[in] key The \e key for the state to access.
  //!
  //! \return Accessor to the state with \e key.
  template <typename T, size_t NS = flecsi_user_space>
  decltype(auto) access_state_(const const_string_t && key)
  {
    return data_t::instance().template dense_accessor<T, NS>( 
      key, mesh_.runtime_id() );
  } // access_state_

  //! \brief Access state registered with type \e T.
  //!
  //! \tparam T All state variables of this type will be returned.
  //! \tparam NS Namespace to use.
  //!
  //! \return A vector of accessors to state registered with type \e T.
  template <typename T, size_t NS = flecsi_user_space>
  decltype(auto) access_type_()
  {
    return data_t::instance().template dense_accessors<T, NS>( 
      mesh_.runtime_id() );
  } // access_type_

  //! \brief Access state registered with type \e T that matches predicate
  //!   function of type \e P.
  //!
  //! \tparam T All state variables of this type will be returned that match
  //!   predicate \e P will be returned.
  //! \tparam P Predicate function type.
  //!
  //! \param[in] predicate Predicate function.
  //!
  //! \return Accessors to the state variables of type \e T matching the
  //!   predicate function.
  template <typename T, typename P>
  decltype(auto) access_type_if_(P && predicate)
  {
    return data_t::instance().template dense_accessors<T, P>(
      std::forward<P>(predicate), mesh_.runtime_id() );
  } // access_type_if


  //============================================================================
  // Global Accessors
  //============================================================================

  //! \brief Register state for the named variable with attributes.
  //!
  //! \tparam T The type of the underlying data to access.
  //!
  //! \param[in] key A name for the state variable, e.g., "density".
  //! \param[in] attributes A bitfield specifying various attributes of the state.
  //!
  //! \return An accessor to the newly registered state.
  template <typename T>
  decltype(auto) register_global_state_(
    const const_string_t && key,
    typename bitfield_t::field_type_t attributes = 0x0 ) 
  {
    return data_t::instance().template register_global_state<T>(
      key, mesh_.runtime_id(), attachment_site_t::global, attributes
    );
  } // register_state_

  //! \brief Access state associated with \e key.
  //!
  //! \tparam T Data type of underlying state.
  //! \tparam NS The namespace in which the state variable is registered.
  //!   See \ref data_t::register_state for additional information.
  //!
  //! \param[in] key The \e key for the state to access.
  //!
  //! \return Accessor to the state with \e key.
  template <typename T, size_t NS = flecsi_user_space>
  decltype(auto) access_global_state_(const const_string_t && key)
  {
    return data_t::instance().template global_accessor<T, NS>(key, mesh_.runtime_id());
  } // access_state_

  //! \brief Access state registered with type \e T.
  //!
  //! \tparam T All state variables of this type will be returned.
  //! \tparam NS Namespace to use.
  //!
  //! \return A vector of accessors to state registered with type \e T.
  template <typename T, size_t NS = flecsi_user_space>
  decltype(auto) access_global_type_()
  {
    return data_t::instance().template global_accessors<T, NS>();
  } // access_type_

  //! \brief Access state registered with type \e T that matches predicate
  //!   function of type \e P.
  //!
  //! \tparam T All state variables of this type will be returned that match
  //!   predicate \e P will be returned.
  //! \tparam P Predicate function type.
  //!
  //! \param[in] predicate Predicate function.
  //!
  //! \return Accessors to the state variables of type \e T matching the
  //!   predicate function.
  template <typename T, typename P>
  decltype(auto) access_global_type_if_(P && predicate)
  {
    return data_t::instance().template global_accessors<T, P>(
      std::forward<P>(predicate));
  } // access_type_if


  //============================================================================
  // General Accessors
  //============================================================================

  //! \brief Return the attributes for the state with \e key.
  //!
  //! \param[in] key The \e key for the state to return attributes for.
  //!
  //! \return The attributes for the state with \e key.
  decltype(auto) state_attributes_(const const_string_t && key)
  {
    return data_t::instance().template meta_data<>((key)).attributes;
  } // state_attribtutes_


  //============================================================================
  // Constructors
  //============================================================================

  //! Default constructor
  burton_mesh_t() {}

  //! Assignment operator (default)
  burton_mesh_t & operator=(const burton_mesh_t &) = default;

  //! Copy constructor
  burton_mesh_t(burton_mesh_t &src);

  //! allow move construction
  burton_mesh_t( burton_mesh_t && ) = default;

  //! move assignment
  burton_mesh_t & operator=(burton_mesh_t && other)
  {
    // store old runtime id
    auto runtime_id = other.mesh_.runtime_id();
    // move the mesh
    mesh_ = std::move( other.mesh_ );
    // move the data
    data_t::instance().move( runtime_id, mesh_.runtime_id() );
    // reset each entity mesh pointer
    for ( auto v : vertices() ) v->reset( mesh_ );
    for ( auto e : edges() ) e->reset( mesh_ );
    for ( auto c : cells() ) c->reset( mesh_ );
    for ( auto c : corners() ) c->reset( mesh_ );
    for ( auto w : wedges() ) w->reset( mesh_ );
    // return mesh
    return *this;
  };

  //! Destructor
  ~burton_mesh_t() 
  {
    // multiple initializations of the data singleton
    data_t::instance().reset( mesh_.runtime_id() );
  }


  //============================================================================
  // Accessors
  //============================================================================

  //! \brief Return the topological dimension of the burton mesh.
  //!
  //! \return A non-negative number describing the highest dimension
  //!   of the entities in the burton mesh, e.g., 3 for a three-dimensional
  //!   burton mesh.
  static constexpr auto num_dimensions()
  {
    return mesh_traits_t::dimension;
  } // dimension
  
  //! \brief Return the time associated with the mesh
  auto time()
  {
    auto soln_time = access_global_state_<real_t, flecsi_internal>( "time" );
    return *soln_time;
  }

  //! \brief Set the time associated with the mesh
  void set_time(real_t soln_time)
  {
    access_global_state_<real_t, flecsi_internal>( "time" ) = soln_time;
  }


  //! \brief Set the time associated with the mesh
  auto increment_time(real_t delta_time)
  {
    auto soln_time = access_global_state_<real_t, flecsi_internal>( "time" );
    (*soln_time) += delta_time;
    return *soln_time;
  }

  //! \brief Return the time associated with the mesh
  auto time_step_counter()
  {
    auto step = access_global_state_<size_t, flecsi_internal>( "time_step" );
    return *step;
  }

  //! \brief Set the time associated with the mesh
  auto increment_time_step_counter(size_t delta = 1)
  {
    auto step = access_global_state_<size_t, flecsi_internal>( "time_step" );
    (*step) += delta;
    return *step;
  }

  //============================================================================
  // Vertex Interface
  //============================================================================

  //! \brief Return number of vertices in the burton mesh.
  //! \return The number of vertices in the burton mesh.
  size_t num_vertices() const
  {
    return mesh_.template num_entities<0, 0>();
  } // num_vertices

  //! \brief Return all vertices in the burton mesh.
  //! \return Return all vertices in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto vertices() const 
  { 
    return mesh_.template entities<0, 0>(); 
  } // vertices

  //! \brief Return vertices associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertices for.
  //!
  //! \param[in] e instance of entity to return vertices for.
  //!
  //! \return Return vertices associated with entity instance \e e as a
  //!    sequence.
  template <class E>
  auto vertices(E * e) const
  {
    return mesh_.template entities<0, 0>(e);
  } // vertices

  //! \brief Return vertices for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get vertices for.
  //!
  //! \param[in] e Entity to get vertices for.
  //!
  //! \return Vertices for entity \e e in domain \e M.
  template <size_t M, class E>
  auto vertices(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<0, M, 0>(e.entity());
  }

  //! \brief Return ids for all vertices in the burton mesh.
  //!
  //! \return Ids for all vertices in the burton mesh.
  auto vertex_ids() const
  {
    return mesh_.template entity_ids<0, 0>();
  } // vertex_ids

  //! \brief Return vertex ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertex ids for.
  //!
  //! \param[in] e instance of entity to return vertex ids for.
  //!
  //! \return Return vertex ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  auto vertex_ids(E * e) const
  {
    return mesh_.template entity_ids<0, 0>(e);
  } // vertex_ids

  //! \brief Return vertices for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get vertices for.
  //!
  //! \param[in] e Entity to get vertices for.
  //!
  //! \return Vertices for entity \e e in domain \e M.
  template <size_t M, class E>
  auto vertex_ids(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entity_ids<0, M, 0>(e.entity());
  }

  //! \brief Return boundary  vertices in the burton mesh.
  //! \return Return all boundary vertices in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto boundary_vertices() const
  { 
    auto verts = vertices();
    auto bnd_verts = verts.filter(
      [](auto v){ return v->is_boundary(); }
    );
    return bnd_verts;
  }


  //============================================================================
  // Edge Interface
  //============================================================================

  //! \brief Return the number of burton mesh edges.
  //! \return The number of burton mesh edges.
  size_t num_edges() const
  {
    return mesh_.template num_entities<1, 0>();
  } // num_edges

  //! \brief Return all edges in the burton mesh.
  //! \return Return all edges in the burton mesh as a sequence for use, e.g., in
  //!   range based for loops.
  auto edges() const { return mesh_.template entities<1, 0>(); } // edges

  //! \brief Return edges for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get edges for.
  //!
  //! \param[in] e Entity to get edges for.
  //!
  //! \return Edges for entity \e e in domain \e M.
  template <size_t M, class E>
  auto edges(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<1, M, 0>(e.entity());
  } // edges

  //! \brief Return ids for all edges in the burton mesh.
  //!
  //! \return Ids for all edges in the burton mesh.
  auto edge_ids() const
  {
    return mesh_.template entity_ids<1, 0>();
  } // edge_ids

  //! \brief Return edge ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return edge ids for.
  //!
  //! \param[in] e instance of entity to return edge ids for.
  //!
  //! \return Return edge ids associated with entity instance \e e as a sequence.
  template <class E>
  auto edge_ids(E * e) const
  {
    return mesh_.template entity_ids<1, 0>(e);
  } // edge_ids

  //! \brief Return boundary edges in the burton mesh.
  //! \return Return all boundary vertices in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto boundary_edges() const
  { 
    auto es = edges();
    auto bnd_es = es.filter(
      [](auto e){ return e->is_boundary(); }
    );
    return bnd_es;
  }

  //============================================================================
  // Cell Interface
  //============================================================================

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  size_t num_cells() const
  {
    return mesh_.template num_entities<num_dimensions(), 0>();
  } // num_cells

  //! \brief Return all cells in the burton mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto cells() const
  {
    return mesh_.template entities<num_dimensions(), 0>();
  } // cells

  //! \brief Return all cells in the burton mesh.
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto cells() // FIXME const
  {
    return mesh_.template entities<num_dimensions(), 0>();
  } // cells

  //! \brief Return cells associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return cells for.
  //!
  //! \param[in] e instance of entity to return cells for.
  //!
  //! \return Return cells associated with entity instance \e e as a sequence.
  template <class E>
  auto cells(E * e) const
  {
    return mesh_.template entities<num_dimensions(), 0>(e);
  } // cells

  //! \brief Return cells for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get cells for.
  //!
  //! \param[in] e Entity to get cells for.
  //!
  //! \return Cells for entity \e e in domain \e M.
  template <size_t M, class E>
  auto cells(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<num_dimensions(), M, 0>(e.entity());
  } // cells

  //! \brief Return ids for all cells in the burton mesh.
  //! \return Ids for all cells in the burton mesh.
  auto cell_ids() const
  {
    return mesh_.template entity_ids<num_dimensions(), 0>();
  } // cell_ids

  //! \brief Return cell ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return cell ids for.
  //!
  //! \param[in] e instance of entity to return cell ids for.
  //!
  //! \return Return cell ids associated with entity instance \e e as a sequence.
  template <class E>
  auto cell_ids(E * e) const
  {
    return mesh_.template entity_ids<num_dimensions(), 0>(e);
  } // cell_ids

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  auto cell_types() const
  {
    auto cs = cells();
    using cell_type_t = decltype( cs[0]->type() );
    std::set< cell_type_t > cell_types;
    for ( auto c : cs ) cell_types.insert( c->type() );
    return cell_types;
  } // num_cells


  //============================================================================
  // Wedge Interface
  //============================================================================

  //! \brief Return number of wedges in the burton mesh.
  //! \return The number of wedges in the burton mesh.
  size_t num_wedges() const
  {
    return mesh_.template num_entities<num_dimensions(), 1>();
  } // num_wedges

  //! \brief Return all wedges in the burton mesh.
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto wedges() // FIXME const
  {
    return mesh_.template entities<num_dimensions(), 1>();
  } // wedges

  //! \brief Return wedges associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return wedges for.
  //!
  //! \param[in] e instance of entity to return wedges for.
  //!
  //! \return Return wedges associated with entity instance \e e as a sequence.
  template <class E>
  auto wedges(E * e) const
  {
    return mesh_.template entities<num_dimensions(), 1>(e);
  } // wedges

  //! \brief Return wedges for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get wedges for.
  //!
  //! \param[in] e Entity to get wedges for.
  //!
  //! \return Wedges for entity \e e in domain \e M.
  template<size_t M, class E>
  auto wedges(flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<num_dimensions(), M, 1>(e.entity());
  }

  //! \brief Return ids for all wedges in the burton mesh.
  //! \return Ids for all wedges in the burton mesh.
  auto wedge_ids() const
  {
    return mesh_.template entity_ids<num_dimensions(), 1>();
  } // wedge_ids

  //! \brief Return wedge ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return wedge ids for.
  //!
  //! \param[in] e instance of entity to return wedge ids for.
  //!
  //! \return Return wedge ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  auto wedge_ids(E * e) const
  {
    return mesh_.template entity_ids<num_dimensions(), 1>(e);
  } // wedge_ids

  //============================================================================
  // Corner Interface
  //============================================================================

  //! \brief Return number of corners in the burton mesh.
  //! \return The number of corners in the burton mesh.
  size_t num_corners() const
  {
    return mesh_.template num_entities<1, 1>();
  } // num_corners

  //! \brief Return all corners in the burton mesh.
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto corners()  // FIXME const
  {
    return mesh_.template entities<1, 1>();
  } // corners

  //! \brief Return corners associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return corners for.
  //!
  //! \param[in] e instance of entity to return corners for.
  //!
  //! \return Return corners associated with entity instance \e e as a sequence.
  template <class E>
  auto corners(E * e) const
  {
    return mesh_.template entities<1, 1>(e);
  } // corners

  //! \brief Return corners for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get corners for.
  //!
  //! \param[in] e Entity to get corners for.
  //!
  //! \return Corners for entity \e e in domain \e M.
  template<size_t M, class E>
  auto corners(flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<1, M, 1>(e.entity());
  }

  //! \brief Return ids for all corners in the burton mesh.
  //! \return Ids for all corners in the burton mesh.
  auto corner_ids() const
  {
    return mesh_.template entity_ids<1, 1>();
  } // corner_ids

  //! \brief Return corner ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return corner ids for.
  //!
  //! \param[in] e instance of entity to return corner ids for.
  //!
  //! \return Return corner ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  auto corner_ids(E * e) const
  {
    return mesh_.template entity_ids<1, 1>(e);
  } // corner_ids


  //============================================================================
  // Region Interface
  //============================================================================

  //! \brief Return the number of regions in the burton mesh.
  //! \return The number of regions in the burton mesh.
  auto num_regions()
  {
    auto n = access_global_state_<size_t, flecsi_internal>( "num_regions" );
    return *n;
  } // num_cells

  //! \brief Return the number of regions in the burton mesh.
  //! \return The number of regions in the burton mesh.
  void set_num_regions(size_t n)
  {
    access_global_state_<size_t, flecsi_internal>( "num_regions" ) = n;
  } // num_cells

  //! \brief Return all cells in the regions mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto regions()
  {
    // get the number of regions
    auto n = num_regions();
    // get the cells
    auto cs = cells();
    // create storage for regions
    using set_type_t = decltype( 
      mesh_.template entities<num_dimensions(), 0>() 
    );

    // create a function to scatter the cells
    using domain_entity_t = decltype( cells()[0] );
    std::function< size_t(domain_entity_t c) > func = 
      [](domain_entity_t c){ return c->region(); };

    // now filter out the cells of each region
    auto region_cells = cs.scatter( func );

    return region_cells;
  } // cells

  //============================================================================
  // Mesh Creation
  //============================================================================

  //! \brief Create a vertex in the burton mesh.
  //!
  //! \param[in] pos The position (coordinates) for the vertex.
  //!
  //! \return Pointer to a vertex created at \e pos.
  vertex_t * create_vertex(const point_t & pos)
  {
    auto p = access_state_<point_t, flecsi_internal>("coordinates");
    p[num_vertices()] = pos;

    auto v = mesh_.template make<vertex_t>( mesh_ );
    mesh_.template add_entity<0, 0>(v);

    return v;
  }


  //! \brief Create a cell in the burton mesh.
  //!
  //! \param[in] verts The vertices defining the cell.
  //!
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell_(V && verts)
  {
    
    cell_t * c;

    switch ( verts.size() ) {
    case (1,2):
      raise_runtime_error( "can't have <3 vertices" );
    case (3):
      c = mesh_.template make< typename mesh_types_t::triangle_cell_t >(mesh_);
      break;
    case (4):
      c = mesh_.template make< typename mesh_types_t::quadrilateral_cell_t >(mesh_);
      break;
    default:
      c = mesh_.template make< typename mesh_types_t::polygonal_cell_t >(mesh_);
      break;
    }

    mesh_.template add_entity<num_dimensions(), 0>(c);

    // FIXME: Need to make mesh interface more general
    mesh_.template init_cell<0>( c, std::forward<V>(verts) );
    return c;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //!
  //! \param[in] verts The vertices defining the cell.
  //!
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell(V && verts)
  {
    return create_cell_( verts );
  } // create_cell

  
  //! \brief Create a cell in the burton mesh.
  //!
  //! \param[in] verts The vertices defining the cell.
  //!
  //! \return Pointer to cell created with \e verts.
  auto  create_cell( std::initializer_list<vertex_t *> verts ) {
    return create_cell_( verts );
  }

  //! \brief Dump the burton mesh to standard out.
  void dump()
  {
    mesh_.dump();
  }
    


  //! \brief Initialize burton mesh state for the number of \e vertices.
  //!
  //! \param[in] vertices The number of \e vertices to initialize the burton mesh
  //!   with.
  void init_parameters(size_t num_nodes)
  {
    // register coordinate state
    data_t::instance().template register_state<point_t, flecsi_internal>(
      "coordinates", num_nodes, mesh_.runtime_id(), 
      attachment_site_t::vertices, persistent
    );
  } // init_parameters

  //! \brief Initialize the burton mesh.
  void init()
  {

    mesh_.template init<0>();
    mesh_.template init_bindings<1>();

    //mesh_.dump();
    
    // Initialize corners
    for (auto cn : corners()) {

      auto cl = cells(cn).front();
      auto es = edges(cn);
      auto vt = vertices(cn).front();

      cn->set_cell(cl);
      cn->add_edge(es.front());
      cn->add_edge(es.back());
      cn->set_vertex(vt);
      
      auto w1 = new wedge_t( mesh_ );
      w1->set_cell(cl);
      w1->set_edge(es.front());
      w1->set_vertex(vt);
      mesh_.template add_entity<num_dimensions(), 1>( w1 );

      auto w2 = new wedge_t( mesh_ );
      w2->set_cell(cl);
      w2->set_edge(es.back());
      w2->set_vertex(vt);
      mesh_.template add_entity<num_dimensions(), 1>( w2 );

      cn->add_wedge(w1);
      cn->add_wedge(w2);
    } // for

    // get the data instance
    data_t & data_ = data_t::instance();

    // register time state
    auto soln_time = 
      data_.template register_global_state<real_t, flecsi_internal>(
        "time", mesh_.runtime_id(), attachment_site_t::global, persistent
      );
    soln_time = 0;

    auto step = data_.template register_global_state<real_t, flecsi_internal>(
      "time_step", mesh_.runtime_id(), attachment_site_t::global, persistent
    );
    *step = 0;

    // register some flags for identifying boundarys
    auto point_flags = data_.template register_state<bitfield_t, flecsi_internal>(
      "point_flags", num_vertices(), mesh_.runtime_id(), 
      attachment_site_t::vertices, persistent
    );
    auto edge_flags = data_.template register_state<bitfield_t, flecsi_internal>(
      "edge_flags", num_edges(), mesh_.runtime_id(), 
      attachment_site_t::edges, persistent
    );

    // now set the boundary flags
    for ( auto e : edges() ) {
      // get the points and cells attached to the edge
      auto points = vertices(e);
      auto zones = cells(e);
      // if there is only one cell, it is a boundary
      if ( zones.size() == 1 ) {
        edge_flags[e].set( 1 << bits::boundary );
        point_flags[ points.front() ].setbit( bits::boundary );
        point_flags[ points.back () ].setbit( bits::boundary );
      }
    } // for

    // identify the cell regions
    auto cell_region = data_.template register_state<size_t, flecsi_internal>(
      "cell_region", num_cells(), mesh_.runtime_id(), 
      attachment_site_t::cells, persistent
    );

    for ( auto c : cells() )
      cell_region[c] = 0;

    // list the regions
    auto num_regions = 
      data_.template register_global_state<size_t, flecsi_internal>(
        "num_regions", mesh_.runtime_id(), 
        attachment_site_t::global, persistent
      );
    *num_regions = 1;

  } // init

  //============================================================================
  // Operators
  //============================================================================

  //! Print some statistics.
  template< std::size_t M >
  friend std::ostream& operator<< (std::ostream& stream, const burton_mesh_t<M>& mesh);

  //============================================================================
  // Private Data 
  //============================================================================

 private:

  mesh_t mesh_;


}; // class burton_mesh_t


////////////////////////////////////////////////////////////////////////////////
// External Class Definitions
////////////////////////////////////////////////////////////////////////////////


//==============================================================================
// Constructors
//==============================================================================

//! \brief Copy constructor.
//! 
//! \param[in] the mesh to copy.
template< std::size_t N >
inline
burton_mesh_t<N>::burton_mesh_t(burton_mesh_t &src)
{
  
  // FIXME!!!
  auto & src_mesh = const_cast<burton_mesh_t&>( src );

  std::vector<vertex_t*> vs;

  init_parameters( src_mesh.num_vertices() );

  // create vertices
  for ( auto v : src_mesh.vertices() ) {
    auto vert = create_vertex( v->coordinates() );
    vs.emplace_back( std::move(vert) );
  }

  // create cells
  for ( auto c : src_mesh.cells() ) {
    auto ids = src_mesh.vertex_ids( c );
    auto n = ids.size();
    std::vector<vertex_t*> elem_vs( n );
    for ( auto i=0; i<n; i++ ) elem_vs[i] = vs[ ids[i] ];
    create_cell( elem_vs );   
  } // for

  // initialize everything
  init();

  auto src_cells = src.cells();

  // override the region ids
  auto num_reg = src.num_regions();
  for ( auto c : cells() ) 
    c->region() = src_cells[c.id()]->region();
  set_num_regions( num_reg );


}

//==============================================================================
// Friends
//==============================================================================

//!  \brief Print some statistics.
//! 
//!  \param[in] stream The stream to print to.
//! 
//!  \param[in] mesh   The mesh object to print stats for.
//! 
//!  \return the stream operator.
template< std::size_t M >
inline
std::ostream& operator<< (std::ostream& stream, const burton_mesh_t<M>& mesh)
{
  using std::endl;
  stream << "Burton mesh:" << endl;
  stream << " + Num Points = " << mesh.num_vertices() << endl;
  stream << " + Num Edges  = " << mesh.num_edges() << endl;
  stream << " + Num Cells  = " << mesh.num_cells() << endl;
  return stream;
}

} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
