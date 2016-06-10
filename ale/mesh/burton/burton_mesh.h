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
#include <sstream>

//! user includes
#include "flecsi/data/data.h"
#include "flecsi/execution/task.h"

#include "ale/mesh/burton/burton_mesh_topology.h"
#include "ale/mesh/burton/burton_types.h"
#include "ale/utils/errors.h"

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

  //! the number of dimensions
  static constexpr auto num_dimensions = mesh_traits_t::num_dimensions;

  //! Type for storing instance of template specialized low level mesh.
  using mesh_topology_t = burton_mesh_topology_t< num_dimensions >;

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
  using face_t = typename mesh_types_t::face_t;

  //! Cell type.
  using cell_t = typename mesh_types_t::cell_t;

  //! Wedge type.
  using wedge_t = typename mesh_types_t::wedge_t;

  //! Corner type.
  using corner_t = typename mesh_types_t::corner_t;

  //! \brief The locations of different bits that we set as flags
  using bits = typename mesh_traits_t::bits;

  //! \brief the type of id for marking boundaries
  using boundary_id_t = typename mesh_traits_t::boundary_id_t;
  using boundary_id_vector_t = typename mesh_traits_t::boundary_id_vector_t;


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
        return 
          data_.template register_state<T>(
            key, num_vertices(), mesh_.runtime_id(), 
            attachment_site_t::vertices, attributes 
          );
        break;
      case attachment_site_t::edges:
        return 
          data_.template register_state<T>(
            key, num_edges(), mesh_.runtime_id(), 
            attachment_site_t::edges, attributes
          );
        break;
      case attachment_site_t::faces:
        return 
          data_.template register_state<T>(
            key, num_faces(), mesh_.runtime_id(), 
            attachment_site_t::faces, attributes 
          );
        break;
      case attachment_site_t::cells:
        return 
          data_.template register_state<T>(
            key, num_cells(), mesh_.runtime_id(), 
            attachment_site_t::cells, attributes 
          );
        break;
      case attachment_site_t::corners:
        return
          data_.template register_state<T>(
            key, num_corners(), mesh_.runtime_id(), 
            attachment_site_t::corners, attributes
          );
        break;
      case attachment_site_t::wedges:
        return 
          data_.template register_state<T>(
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
    for ( auto f : faces() ) f->reset( mesh_ );
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
    return mesh_.template num_entities<vertex_t::dimension, vertex_t::domain>();
  } // num_vertices

  //! \brief Return all vertices in the burton mesh.
  //! \return Return all vertices in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto vertices() const 
  { 
    return mesh_.template entities<vertex_t::dimension, vertex_t::domain>(); 
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
    return mesh_.template entities<vertex_t::dimension, vertex_t::domain>(e);
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
    return mesh_.template entities<vertex_t::dimension, M, vertex_t::domain>(e.entity());
  }


  //! \brief Return vertices associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertices for.
  //!
  //! \param[in] e instance of entity to return vertices for.
  //!
  //! \return Return vertices associated with entity instance \e e as a
  //!    sequence.
  template <
    typename P,
    typename = typename std::enable_if_t< utils::is_callable_v<P> >
  >
  auto vertices( P && p ) const
  {
    
    auto vs = vertices();
    return vs.filter( std::forward<P>(p) );
  } // vertices

  //! \brief Return ids for all vertices in the burton mesh.
  //!
  //! \return Ids for all vertices in the burton mesh.
  auto vertex_ids() const
  {
    return mesh_.template entity_ids<vertex_t::dimension, vertex_t::domain>();
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
    return mesh_.template entity_ids<vertex_t::dimension, vertex_t::domain>(e);
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
    return mesh_.template entity_ids<vertex_t::dimension, M, vertex_t::domain>(e.entity());
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
    return mesh_.template num_entities<edge_t::dimension, edge_t::domain>();
  } // num_edges

  //! \brief Return all edges in the burton mesh.
  //! \return Return all edges in the burton mesh as a sequence for use, e.g., in
  //!   range based for loops.
  auto edges() const { 
    return mesh_.template entities<edge_t::dimension, 0>(); 
  } // edges

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
    return mesh_.template entities<edge_t::dimension, M, edge_t::domain>(e.entity());
  } // edges

  //! \brief Return ids for all edges in the burton mesh.
  //!
  //! \return Ids for all edges in the burton mesh.
  auto edge_ids() const
  {
    return mesh_.template entity_ids<edge_t::dimension, edge_t::domain>();
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
    return mesh_.template entity_ids<edge_t::dimension, edge_t::domain>(e);
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
  // Face Interface
  //============================================================================

  //! \brief Return the number of faces in the burton mesh.
  //! \return The number of faces in the burton mesh.
  size_t num_faces() const
  {
    return mesh_.template num_entities<face_t::dimension, face_t::domain>();
  } // num_faces

  //! \brief Return all faces in the burton mesh.
  //!
  //! \return Return all faces in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto faces() const
  {
    return mesh_.template entities<face_t::dimension, face_t::domain>();
  } // faces

  //! \brief Return all faces in the burton mesh.
  //! \return Return all faces in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto faces() // FIXME const
  {
    return mesh_.template entities<face_t::dimension, face_t::domain>();
  } // faces

  //! \brief Return faces associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return faces for.
  //!
  //! \param[in] e instance of entity to return faces for.
  //!
  //! \return Return faces associated with entity instance \e e as a sequence.
  template <class E>
  auto faces(E * e) const
  {
    return mesh_.template entities<face_t::dimension, face_t::domain>(e);
  } // faces

  //! \brief Return faces for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get faces for.
  //!
  //! \param[in] e Entity to get faces for.
  //!
  //! \return Faces for entity \e e in domain \e M.
  template <size_t M, class E>
  auto faces(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<face_t::dimension, M, face_t::domain>(e.entity());
  } // faces

  //! \brief Return ids for all faces in the burton mesh.
  //! \return Ids for all faces in the burton mesh.
  auto face_ids() const
  {
    return mesh_.template entity_ids<face_t::dimension, face_t::domain>();
  } // face_ids

  //! \brief Return face ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return face ids for.
  //!
  //! \param[in] e instance of entity to return face ids for.
  //!
  //! \return Return face ids associated with entity instance \e e as a sequence.
  template <class E>
  auto face_ids(E * e) const
  {
    return mesh_.template entity_ids<face_t::dimension, face_t::domain>(e);
  } // face_ids

  //============================================================================
  // Cell Interface
  //============================================================================

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  size_t num_cells() const
  {
    return mesh_.template num_entities<cell_t::dimension, cell_t::domain>();
  } // num_cells

  //! \brief Return all cells in the burton mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto cells() const
  {
    return mesh_.template entities<cell_t::dimension, cell_t::domain>();
  } // cells

  //! \brief Return all cells in the burton mesh.
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto cells() // FIXME const
  {
    return mesh_.template entities<cell_t::dimension, cell_t::domain>();
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
    return mesh_.template entities<cell_t::dimension, cell_t::domain>(e);
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
    return mesh_.template entities<cell_t::dimension, M, cell_t::domain>(e.entity());
  } // cells

  //! \brief Return ids for all cells in the burton mesh.
  //! \return Ids for all cells in the burton mesh.
  auto cell_ids() const
  {
    return mesh_.template entity_ids<cell_t::dimension, cell_t::domain>();
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
    return mesh_.template entity_ids<cell_t::dimension, cell_t::domain>(e);
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
    return mesh_.template num_entities<wedge_t::dimension, wedge_t::domain>();
  } // num_wedges

  //! \brief Return all wedges in the burton mesh.
  //!
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto wedges() const
  {
    return mesh_.template entities<wedge_t::dimension, wedge_t::domain>();
  } // wedges

  //! \brief Return all wedges in the burton mesh.
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto wedges() // FIXME const
  {
    return mesh_.template entities<wedge_t::dimension, wedge_t::domain>();
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
    return mesh_.template entities<wedge_t::dimension, wedge_t::domain>(e);
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
  auto wedges(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<wedge_t::dimension, M, wedge_t::domain>(e.entity());
  }

  //! \brief Return ids for all wedges in the burton mesh.
  //! \return Ids for all wedges in the burton mesh.
  auto wedge_ids() const
  {
    return mesh_.template entity_ids<wedge_t::dimension, wedge_t::domain>();
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
    return mesh_.template entity_ids<wedge_t::dimension, wedge_t::domain>(e);
  } // wedge_ids

  //============================================================================
  // Corner Interface
  //============================================================================

  //! \brief Return number of corners in the burton mesh.
  //! \return The number of corners in the burton mesh.
  size_t num_corners() const
  {
    return mesh_.template num_entities<corner_t::dimension, corner_t::domain>();
  } // num_corners

  //! \brief Return all corners in the burton mesh.
  //!
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto corners() const
  {
    return mesh_.template entities<corner_t::dimension, corner_t::domain>();
  } // corners

  //! \brief Return all corners in the burton mesh.
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  auto corners() // FIXME const
  {
    return mesh_.template entities<corner_t::dimension, corner_t::domain>();
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
    return mesh_.template entities<corner_t::dimension, corner_t::domain>(e);
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
  auto corners(const flecsi::domain_entity<M, E> & e) const
  {
    return mesh_.template entities<corner_t::dimension, M, corner_t::domain>(e.entity());
  }

  //! \brief Return ids for all corners in the burton mesh.
  //! \return Ids for all corners in the burton mesh.
  auto corner_ids() const
  {
    return mesh_.template entity_ids<corner_t::dimension, corner_t::domain>();
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
    return mesh_.template entity_ids<corner_t::dimension, corner_t::domain>(e);
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

  //! \brief set the number of regions in the burton mesh.
  //! \param [in]  n  The number of regions in the burton mesh.
  void set_num_regions(size_t n)
  {
    access_global_state_<size_t, flecsi_internal>( "num_regions" ) = n;
  } // num_cells


  //! \brief Return the number of regions in the burton mesh.
  //! \param [in]  n  The number of regions in the burton mesh.
  template< typename T >
  void set_regions(T * region_ids)
  {
    for ( auto c : cells() )
      c->region() = region_ids[c.id()];
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
      mesh_.template entities<cell_t::dimension, cell_t::domain>() 
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
  // Element Creation
  //============================================================================


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell(
    V && verts,
    typename std::enable_if_t< 
      std::is_same_v< typename std::decay_t<V>::value_type, vertex_t* > &&
      std::remove_pointer_t<typename std::decay_t<V>::value_type>::num_dimensions == 2
    >* = nullptr ) 
  {
    return create_2d_element_from_verts_<cell_t>( std::forward<V>(verts) );
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto  create_cell( 
    std::initializer_list<V*> verts,
    typename std::enable_if_t< 
      std::is_same_v<V, vertex_t> && V::num_dimensions == 2 
    >* = nullptr ) 
  {
    return create_2d_element_from_verts_<cell_t>( verts );
  }


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell(
    V && verts,
    typename std::enable_if_t< 
      std::is_same_v< typename std::decay_t<V>::value_type, vertex_t* > &&
      std::remove_pointer_t<typename std::decay_t<V>::value_type>::num_dimensions == 3
    >* = nullptr ) 
  {
    return create_3d_element_from_verts_( std::forward<V>(verts) );
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto  create_cell( 
    std::initializer_list<V*> verts,
    typename std::enable_if_t< 
      std::is_same_v<V, vertex_t> && V::num_dimensions == 3
    >* = nullptr ) 
  {
    return create_3d_element_from_verts_( verts );
  }

  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  template< typename F >
  auto create_cell(
    F && faces,
    typename std::enable_if_t< 
      std::is_same_v< typename std::decay_t<F>::value_type, face_t* >  &&
      std::remove_pointer_t<typename std::decay_t<F>::value_type>::num_dimensions == 3
    >* = nullptr ) 
  {
    return create_3d_element_from_faces_( std::forward<F>(faces) );
  } // create_cell
  
  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  auto  create_cell( std::initializer_list<face_t *> faces ) {
    return create_3d_element_from_faces_( faces );
  }



  //! \brief Create a face in the burton mesh.
  //! \param[in] verts The vertices defining the face.
  //! \return Pointer to cell created with \e faces.
  template< typename V >
  auto create_face(V && verts)
  {
    return create_2d_element_from_verts_<face_t>( std::forward<V>(verts) );
  } // create_cell

  
  //! \brief Create a face in the burton mesh.
  //! \param[in] verts The vertices defining the face.
  //! \return Pointer to cell created with \e faces.
  auto  create_face( std::initializer_list<vertex_t *> verts ) {
    return create_2d_element_from_verts_<face_t>( verts );
  }


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
    mesh_.template add_entity<vertex_t::dimension, vertex_t::domain>(v);

    return v;
  }


  //============================================================================
  // Mesh Creation
  //============================================================================

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

  //!---------------------------------------------------------------------------
  //! \brief Initialize the burton mesh.
  //!---------------------------------------------------------------------------
  void init()
  {

    mesh_.template init<0>();
    mesh_.template init_bindings<1>();

    //mesh_.dump();

#if 0
    // make sure faces point from first to second cell
    for(auto f : faces()) {
      auto n = f->normal();
      auto fx = f->centroid();
      auto c = cells(f).front();
      auto cx = c->centroid();
      auto delta = fx - cx;
      auto dot = dot_product( n, delta );
      if ( dot < 0 ) {
        std::cout << "reversing" << std::endl;
        mesh_.template reverse_entities<vertex_t::dimension, vertex_t::domain>(f);
      }
    } // for
#endif

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

    // register some flags for identifying boundarys and various other things
    auto point_flags = data_.template register_state<bitfield_t, flecsi_internal>(
      "point_flags", num_vertices(), mesh_.runtime_id(), 
      attachment_site_t::vertices, persistent
    );
    auto edge_flags = data_.template register_state<bitfield_t, flecsi_internal>(
      "edge_flags", num_edges(), mesh_.runtime_id(), 
      attachment_site_t::edges, persistent
    );

    // register some flags for associating boundaries with entities
    data_.template register_state<boundary_id_vector_t, flecsi_internal>(
      "point_boundary_ids", num_vertices(), mesh_.runtime_id(), 
      attachment_site_t::vertices, persistent
    );
    data_.template register_state<boundary_id_vector_t, flecsi_internal>(
      "edge_boundary_ids", num_edges(), mesh_.runtime_id(), 
      attachment_site_t::edges, persistent
    );
    data_.template register_state<boundary_id_t, flecsi_internal>(
      "face_boundary_ids", num_faces(), mesh_.runtime_id(), 
      attachment_site_t::faces, persistent
    );

    // now set the boundary flags.
    for ( auto f : faces() ) {
      // get the points and cells attached to the edge
      auto ps = vertices(f);
      auto es = edges(f);
      auto cs = cells(f);
      // if there is only one cell, it is a boundary
      if ( f->is_boundary() ) {
        for ( auto e : es ) 
          edge_flags[e].setbit( bits::boundary );
        for ( auto p : ps ) 
          point_flags[ p ].setbit( bits::boundary );
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




  //!---------------------------------------------------------------------------
  //! \brief Check the burton mesh.
  //!---------------------------------------------------------------------------
  bool is_valid( bool raise_on_error = true )
  {
    // some includes
    using math::dot_product;

    // a lambda function for raising errors or returning 
    // false
    auto raise_or_return = [=]( std::ostream & msg )
      {
        if ( raise_on_error )
          raise_runtime_error( msg.rdbuf() );
        else 
          std::cerr << msg.rdbuf() << std::endl;
        return false;
      };
    
    // we use a stringstream to construct messages and pass them 
    // simultanesously
    std::stringstream ss;


    // make sure face normal points out from first cell
    for(auto f : faces()) {
      auto n = f->normal();
      auto fx = f->midpoint();
      auto c = cells(f).front();
      auto cx = c->midpoint();
      auto delta = fx - cx;
      auto dot = dot_product( n, delta );      
      if ( dot < 0 ) {
        return raise_or_return( ss << "Face " << f.id() << " has opposite normal" );
      }
    } 

    // check all the corners and wedges
    for(auto cn : corners()) {

      auto cs = cells(cn);
      auto fs = faces(cn);
      auto es = edges(cn);
      auto vs = vertices(cn);
      auto ws = wedges(cn);

      if ( cs.size() != 1 ) 
        return raise_or_return( 
          ss << "Corner " << cn.id() << " has " << cs.size() << "/=1 cells" );

      if ( fs.size() != num_dimensions ) 
        return raise_or_return( 
          ss << "Corner " << cn.id() << " has " << fs.size() << "/=" 
          << num_dimensions << " faces" );

      if ( es.size() != num_dimensions ) 
        return raise_or_return( 
          ss << "Corner " << cn.id() << " has " << es.size() << "/=" 
          << num_dimensions << " edges" );

      if ( vs.size() != 1 )
        return raise_or_return( 
          ss << "Corner " << cn.id() << " has " << vs.size() << "/=1 vertices" );

      auto cl = cs.front();
      auto vt = vs.front();
      
      if ( ws.size() % 2 != 0 )
        return raise_or_return( 
          ss << "Corner " << cn.id() << " has " << ws.size() << "%2/=0 wedges" );

      for ( auto wg = ws.begin(); wg != ws.end();  ) 
        for ( auto i=0; i<2 && wg != ws.end(); i++, ++wg)
        {
          auto cls = cells( *wg );
          auto fs = faces( *wg );
          auto es = edges( *wg );
          auto vs = vertices( *wg );
          auto cns = corners( *wg );
          if ( cls.size() != 1 )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has " << cls.size() << "/=1 cells" );
          if ( fs.size() != 1 )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has " << fs.size() << "/=1 faces" );
          if ( es.size() != 1 )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has " << es.size() << "/=1 edges" );
          if ( vs.size() != 1 )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has " << vs.size() << "/=1 vertices" );
          if ( cns.size() != 1 )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has " << cns.size() << "/=1 corners" );          
          auto vert = vs.front();
          auto cell = cls.front();
          auto corn = cns.front();
          if ( vert != vt )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has incorrect vertex " 
              << vert.id() << "!=" << vt.id() );
          if ( cell != cl )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has incorrect cell " 
              << cell.id() << "!=" << cl.id() );
          if ( corn != cn )
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has incorrect corner " 
              << corn.id() << "!=" << cn.id() );
          auto fc = fs.front();            
          auto fx = fc->midpoint();
          auto cx = cl->midpoint();
          auto delta = fx - cx;
          real_t dot;
          if ( i == 0 ) {
            auto n = wg->facet_normal_right();
            dot = dot_product( n, delta );
          }
          else {
            auto n = wg->facet_normal_left();
            dot = dot_product( n, delta );
          }
          if ( dot < 0 ) 
            return raise_or_return( 
              ss << "Wedge " << (*wg).id() << " has opposite normal" );
        } // wedges
      
    } // corners

    return true;

  }


  //============================================================================
  // Boundary conditions
  //============================================================================
  template< typename P >
  boundary_id_t install_boundary( P && p ) 
  {
    // increment the boundary face storage
    auto this_bnd = face_sets_.size();
    auto num_bnd = this_bnd + 1;
    face_sets_.resize( num_bnd );
    edge_sets_.resize( num_bnd );
    vert_sets_.resize( num_bnd );

    auto & this_bnd_faces = face_sets_[ this_bnd ];
    auto & this_bnd_edges = edge_sets_[ this_bnd ];
    auto & this_bnd_verts = vert_sets_[ this_bnd ];

    // add the face tags and collect the attached edge and vertices
    for ( auto f : faces() )
      if ( p( f ) ) {
        auto es = edges( f );
        auto vs = vertices( f );
        f->tag_boundary( this_bnd );
        this_bnd_faces.emplace_back( f );
        this_bnd_edges.reserve( this_bnd_edges.size() + es.size() );
        this_bnd_verts.reserve( this_bnd_verts.size() + vs.size() );
        for ( auto e : es ) 
          this_bnd_edges.emplace_back( e );
        for ( auto v : vs )
          this_bnd_verts.emplace_back( v );
      }

    // need to remove duplicates from edge and vertex lists
    std::sort( this_bnd_edges.begin(), this_bnd_edges.end() );
    std::sort( this_bnd_verts.begin(), this_bnd_verts.end() );

    this_bnd_edges.erase( 
      std::unique( this_bnd_edges.begin(), this_bnd_edges.end() ), 
      this_bnd_edges.end() 
    );
    this_bnd_verts.erase( 
      std::unique( this_bnd_verts.begin(), this_bnd_verts.end() ), 
      this_bnd_verts.end() 
    );

    // add the edge tags
    for ( auto e : this_bnd_edges )
      e->tag_boundary( this_bnd );
    // add the vertex tags
    for ( auto v : this_bnd_verts )
      v->tag_boundary( this_bnd );

    return num_bnd;
  }

  //============================================================================
  // Operators
  //============================================================================

  //! Print some statistics.
  template< std::size_t M >
  friend std::ostream& operator<< (std::ostream& stream, const burton_mesh_t<M>& mesh);



  //============================================================================
  // Private Members
  //============================================================================

 private:


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename E, typename V >
  auto create_2d_element_from_verts_( V && verts  )
  {
    
    E * e;

    switch ( verts.size() ) {
    case (1,2):
      raise_runtime_error( "can't have <3 vertices" );
    case (3):
      e = mesh_.template make< burton_triangle_t<num_dimensions> >(mesh_);
      break;
    case (4):
      e = mesh_.template make< burton_quadrilateral_t<num_dimensions> >(mesh_);
      break;
    default:
      e = mesh_.template make< burton_polygon_t<num_dimensions> >(mesh_);
      break;
    }

    mesh_.template add_entity<E::dimension, E::domain>( e );
    mesh_.template init_entity<E::domain, E::dimension, vertex_t::dimension>( e, std::forward<V>(verts) );
    return e;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_3d_element_from_verts_( V && verts )
  {
    
    cell_t * c;

    switch ( verts.size() ) {
    case (1,2,3):
      raise_runtime_error( "can't have <4 vertices" );
    case (4):
      c = mesh_.template make< burton_tetrahedron_t >(mesh_);
      break;
    case (8):
      c = mesh_.template make< burton_hexahedron_t >(mesh_);
      break;
    default:
      raise_runtime_error( "can't build polyhedron from vertices alone" );      
      break;
    }

    mesh_.template add_entity<cell_t::dimension, cell_t::domain>(c);
    mesh_.template init_cell<cell_t::domain>( c, std::forward<V>(verts) );
    return c;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  template< typename F >
  auto create_3d_element_from_faces_( F && faces )
  {
    
    cell_t * c;

    switch ( faces.size() ) {
    case (1,2,3):
      raise_runtime_error( "can't have <4 vertices" );
    default:
      c = mesh_.template make< burton_polyhedron_t >( mesh_, std::forward<F>(faces) );
      break;
    }

    mesh_.template add_entity<cell_t::dimension, cell_t::domain>(c);
    mesh_.template init_entity<cell_t::domain, cell_t::dimension, face_t::dimension>( c, std::forward<F>(faces) );
    return c;
  } // create_cell


  //============================================================================
  // Private Data 
  //============================================================================

  mesh_topology_t mesh_;

  std::vector< std::vector<face_t*> >   face_sets_;
  std::vector< std::vector<edge_t*> >   edge_sets_;
  std::vector< std::vector<vertex_t*> > vert_sets_;


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
  stream << " + Num Faces  = " << mesh.num_faces() << endl;
  stream << " + Num Cells  = " << mesh.num_cells() << endl;
  return stream;
}


} // namespace mesh
} // namespace ale


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
