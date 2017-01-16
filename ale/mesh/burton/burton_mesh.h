/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "ale/mesh/burton/burton_mesh_topology.h"
#include "ale/mesh/burton/burton_types.h"
#include "ale/utils/errors.h"

#include "flecsi/data/data.h"
#include "flecsi/execution/task.h"

// system includes
#include <set>
#include <string>
#include <sstream>


namespace ale {
namespace mesh {

//! This namespace is used to expose enumerations and types.
namespace burton {
  
////////////////////////////////////////////////////////////////////////////////
/// \brief The burton mesh index spaces.
////////////////////////////////////////////////////////////////////////////////
enum index_spaces_t : size_t {
  vertices,
  edges,
  faces,
  cells,
  corners,
  wedges,
};

////////////////////////////////////////////////////////////////////////////////
/// \brief Attributes for flecsi.
////////////////////////////////////////////////////////////////////////////////
enum data_attributes_t : size_t {
  persistent
};


} // namespace burton

////////////////////////////////////////////////////////////////////////////////
/// \brief A specialization of the flecsi low-level mesh topology, state and
///   execution models.
/// \tparam N The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_mesh_t  : public burton_mesh_topology_t<N>

{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the base type
  using base_t = burton_mesh_topology_t<N>;

  //! \brief the mesh types
  using types_t = burton_types_t<N>;

  //! \brief the mesh traits
  using config_t = typename types_t::config_t;

  //! \brief the number of dimensions
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! a compile string type
  using const_string_t = typename config_t::const_string_t;
  
  //! a bitfield type
  using bitfield_t = typename config_t::bitfield_t;

  //! Integer data type.
  using integer_t = typename config_t::integer_t;

  //! Floating point data type.
  using real_t = typename config_t::real_t;

  //! The size type.
  using size_t = typename config_t::size_t;

  //! The type used for loop indexing
  using counter_t = typename config_t::counter_t;

  //! Point data type.
  using point_t = typename config_t::point_t;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //! Vertex type.
  using vertex_t = typename types_t::vertex_t;

  //! Edge type.
  using edge_t = typename types_t::edge_t;

  //! Cell type.
  using face_t = typename types_t::face_t;

  //! Cell type.
  using cell_t = typename types_t::cell_t;

  //! Wedge type.
  using wedge_t = typename types_t::wedge_t;

  //! Corner type.
  using corner_t = typename types_t::corner_t;

  //! \brief The locations of different bits that we set as flags
  using bits = typename config_t::bits;

  //! \brief the type of id for marking boundaries
  using tag_t = typename config_t::tag_t;
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Default constructor
  burton_mesh_t() = default;

  //! \brief Assignment operator (default)
  burton_mesh_t & operator=(const burton_mesh_t &) = default;

  //! \brief Copy constructor
  burton_mesh_t(const burton_mesh_t &src) {

    std::vector<vertex_t*> vs;

    init_parameters( src.num_vertices() );

    // create vertices
    for ( auto v : src.vertices() ) {
      auto vert = create_vertex( v->coordinates() );
      vs.emplace_back( std::move(vert) );
    }

    // create cells
    for ( auto c : src.cells() ) {
      auto verts = src.vertices( c );
      auto n = verts.size();
      std::vector<vertex_t*> elem_vs( n );
      for ( auto i=0; i<n; i++ ) elem_vs[i] = vs[ verts[i].id() ];
        create_cell( elem_vs );   
    } // for

    // initialize everything
    init();

    // get the source cells
    auto src_cells = src.cells();

    // override the region ids
    auto num_reg = src.num_regions();
    for ( auto c : cells() ) 
      c->region() = src_cells[c.id()]->region();
    set_num_regions( num_reg );
  }

  //! \brief allow move construction
  burton_mesh_t( burton_mesh_t && ) = default;

  //! \brief move assignment
  burton_mesh_t & operator=(burton_mesh_t && other)
  {
    // call the base type operator to move the data
    base_t::operator=(std::move(other));
    // reset each entity mesh pointer
    for ( auto v : vertices() ) v->reset( *this );
    for ( auto e : edges() ) e->reset( *this );
    for ( auto f : faces() ) f->reset( *this );
    for ( auto c : cells() ) c->reset( *this );
    for ( auto c : corners() ) c->reset( *this );
    for ( auto w : wedges() ) w->reset( *this );
    // return mesh
    return *this;
  };

  //! Destructor
  ~burton_mesh_t() = default;

  //============================================================================
  //! \brief Return the size of an associated index space/
  //! \param [in] index_space_id
  //============================================================================
  size_t indices( size_t index_space_id ) const override
  {
    switch(index_space_id) {
      case burton::vertices:
        return base_t::num_entities(vertex_t::dimension);
      case burton::edges:
        return base_t::num_entities(edge_t::dimension);
      case burton::faces:
        return base_t::num_entities(face_t::dimension);
      case burton::cells:
        return base_t::num_entities(cell_t::dimension);
      case burton::corners:
        return 
          base_t::template 
            num_entities<corner_t::dimension, corner_t::domain>();
      case burton::wedges:
        return 
          base_t::template num_entities<wedge_t::dimension, wedge_t::domain>();
      default:
        raise_runtime_error("unknown index space");
        return 0;
    } // switch
  } // indices


  //============================================================================
  // Accessors
  //============================================================================

  //! \brief Return the time associated with the mesh
  auto time()
  {
    auto soln_time = get_accessor(*this, mesh, time, real_t, global, 0 );
    return *soln_time;
  }

  //! \brief Set the time associated with the mesh
  //! \param [in] soln_time  The solution time.
  void set_time(real_t soln_time)
  {
    get_accessor(*this, mesh, time, real_t, global, 0 ) = soln_time;
  }


  //! \brief Set the time associated with the mesh
  //! \param [in] delta_time  The solution time increment.
  //! \return The new solution time.
  auto increment_time(real_t delta_time)
  {
    auto soln_time = get_accessor(*this, mesh, time, real_t, global, 0 );
    (*soln_time) += delta_time;
    return *soln_time;
  }

  //! \brief Return the time associated with the mesh
  auto time_step_counter()
  {
    auto step = get_accessor(*this, mesh, time_step, size_t, global, 0 );
    return *step;
  }

  //! \brief Increment the time step counter associated with the mesh
  //! \param [in] delta  The counter increment.
  //! \return The new counter value.
  auto increment_time_step_counter(size_t delta = 1)
  {
    auto step = get_accessor(*this, mesh, time_step, size_t, global, 0 );
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
    return 
      base_t::template num_entities<vertex_t::dimension, vertex_t::domain>();
  }

  //! \brief Return all vertices in the burton mesh.
  //! \return Return all vertices in the burton mesh as a sequence for use,
  //!   e.g., in range based for loops.
  decltype(auto) vertices() const 
  { 
    return base_t::template entities<vertex_t::dimension, vertex_t::domain>(); 
  }

  //! \brief Return vertices associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertices for.
  //!
  //! \param[in] e instance of entity to return vertices for.
  //!
  //! \return Return vertices associated with entity instance \e e as a
  //!    sequence.
  template <class E>
  decltype(auto) vertices(E * e) const
  {
    return base_t::template entities<vertex_t::dimension, vertex_t::domain>(e);
  }

  //! \brief Return vertices for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get vertices for.
  //!
  //! \param[in] e Entity to get vertices for.
  //!
  //! \return Vertices for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) vertices(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return 
      base_t::template entities<vertex_t::dimension, M, vertex_t::domain>(
        e.entity()
      );
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
  decltype(auto) vertices( P && p ) const
  {
    
    auto vs = vertices();
    return vs.filter( std::forward<P>(p) );
  }

  //! \brief Return ids for all vertices in the burton mesh.
  //!
  //! \return Ids for all vertices in the burton mesh.
  decltype(auto) vertex_ids() const
  {
    return base_t::template entity_ids<vertex_t::dimension, vertex_t::domain>();
  }

  //! \brief Return vertex ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertex ids for.
  //!
  //! \param[in] e instance of entity to return vertex ids for.
  //!
  //! \return Return vertex ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) vertex_ids(E * e) const
  {
    return 
      base_t::template entity_ids<vertex_t::dimension, vertex_t::domain>(e);
  }

  //! \brief Return vertices for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get vertices for.
  //!
  //! \param[in] e Entity to get vertices for.
  //!
  //! \return Vertices for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) 
  vertex_ids(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return 
      base_t::template entity_ids<vertex_t::dimension, M, vertex_t::domain>(
        e.entity()
      );
  }

  //! \brief Return boundary  vertices in the burton mesh.
  //! \return Return all boundary vertices in the burton mesh as a sequence for 
  //!         use, e.g., in range based for loops.
  decltype(auto) boundary_vertices() const
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
    return base_t::template num_entities<edge_t::dimension, edge_t::domain>();
  }

  //! \brief Return all edges in the burton mesh.
  //! \return Return all edges in the burton mesh as a sequence for use, e.g., 
  //!         in range based for loops.
  decltype(auto) edges() const { 
    return base_t::template entities<edge_t::dimension, 0>(); 
  }

  //! \brief Return edges for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get edges for.
  //!
  //! \param[in] e Entity to get edges for.
  //!
  //! \return Edges for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) edges(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return base_t::template entities<edge_t::dimension, M, edge_t::domain>(
      e.entity()
    );
  }

  //! \brief Return ids for all edges in the burton mesh.
  //!
  //! \return Ids for all edges in the burton mesh.
  decltype(auto) edge_ids() const
  {
    return base_t::template entity_ids<edge_t::dimension, edge_t::domain>();
  }

  //! \brief Return edge ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return edge ids for.
  //!
  //! \param[in] e instance of entity to return edge ids for.
  //!
  //! \return Return edge ids associated with entity instance \e e as a 
  //!         sequence.
  template <class E>
  decltype(auto) edge_ids(E * e) const
  {
    return base_t::template entity_ids<edge_t::dimension, edge_t::domain>(e);
  }

  //! \brief Return boundary edges in the burton mesh.
  //! \return Return all boundary vertices in the burton mesh as a sequence for 
  //!         use, e.g., in range based for loops.
  decltype(auto) boundary_edges() const
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
    return base_t::template num_entities<face_t::dimension, face_t::domain>();
  } // num_faces

  //! \brief Return all faces in the burton mesh.
  //!
  //! \return Return all faces in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) faces() const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>();
  }

  //! \brief Return all faces in the burton mesh.
  //! \return Return all faces in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) faces() // FIXME const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>();
  }

  //! \brief Return faces associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return faces for.
  //!
  //! \param[in] e instance of entity to return faces for.
  //!
  //! \return Return faces associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) faces(E * e) const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>(e);
  }

  //! \brief Return faces for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get faces for.
  //!
  //! \param[in] e Entity to get faces for.
  //!
  //! \return Faces for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) faces(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return 
      base_t::template entities<face_t::dimension, M, face_t::domain>(
        e.entity()
      );
  }

  //! \brief Return ids for all faces in the burton mesh.
  //! \return Ids for all faces in the burton mesh.
  decltype(auto) face_ids() const
  {
    return base_t::template entity_ids<face_t::dimension, face_t::domain>();
  }

  //! \brief Return face ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return face ids for.
  //!
  //! \param[in] e instance of entity to return face ids for.
  //!
  //! \return Return face ids associated with entity instance \e e as a 
  //!         sequence.
  template <class E>
  decltype(auto) face_ids(E * e) const
  {
    return base_t::template entity_ids<face_t::dimension, face_t::domain>(e);
  }

  //============================================================================
  // Cell Interface
  //============================================================================

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  size_t num_cells() const
  {
    return base_t::template num_entities<cell_t::dimension, cell_t::domain>();
  }

  //! \brief Return all cells in the burton mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) cells() const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>();
  }

  //! \brief Return all cells in the burton mesh.
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) cells() // FIXME const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>();
  }

  //! \brief Return cells associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return cells for.
  //!
  //! \param[in] e instance of entity to return cells for.
  //!
  //! \return Return cells associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) cells(E * e) const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>(e);
  }

  //! \brief Return cells for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get cells for.
  //!
  //! \param[in] e Entity to get cells for.
  //!
  //! \return Cells for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) cells(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return 
      base_t::template entities<cell_t::dimension, M, cell_t::domain>(
        e.entity()
      );
  }

  //! \brief Return ids for all cells in the burton mesh.
  //! \return Ids for all cells in the burton mesh.
  decltype(auto) cell_ids() const
  {
    return base_t::template entity_ids<cell_t::dimension, cell_t::domain>();
  }

  //! \brief Return cell ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return cell ids for.
  //!
  //! \param[in] e instance of entity to return cell ids for.
  //!
  //! \return Return cell ids associated with entity instance \e e as a 
  //!         sequence.
  template <class E>
  decltype(auto) cell_ids(E * e) const
  {
    return base_t::template entity_ids<cell_t::dimension, cell_t::domain>(e);
  }

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  decltype(auto) cell_types() const
  {
    auto cs = cells();
    using cell_type_t = decltype( cs[0]->type() );
    std::set< cell_type_t > cell_types;
    for ( auto c : cs ) cell_types.insert( c->type() );
    return cell_types;
  }


  //============================================================================
  // Wedge Interface
  //============================================================================

  //! \brief Return number of wedges in the burton mesh.
  //! \return The number of wedges in the burton mesh.
  size_t num_wedges() const
  {
    return base_t::template num_entities<wedge_t::dimension, wedge_t::domain>();
  }

  //! \brief Return all wedges in the burton mesh.
  //!
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) wedges() const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>();
  }

  //! \brief Return all wedges in the burton mesh.
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) wedges() // FIXME const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>();
  }

  //! \brief Return wedges associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return wedges for.
  //!
  //! \param[in] e instance of entity to return wedges for.
  //!
  //! \return Return wedges associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) wedges(E * e) const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>(e);
  }

  //! \brief Return wedges for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get wedges for.
  //!
  //! \param[in] e Entity to get wedges for.
  //!
  //! \return Wedges for entity \e e in domain \e M.
  template<size_t M, class E>
  decltype(auto) wedges(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return 
      base_t::template entities<wedge_t::dimension, M, wedge_t::domain>(
        e.entity()
      );
  }

  //! \brief Return ids for all wedges in the burton mesh.
  //! \return Ids for all wedges in the burton mesh.
  decltype(auto) wedge_ids() const
  {
    return base_t::template entity_ids<wedge_t::dimension, wedge_t::domain>();
  }

  //! \brief Return wedge ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return wedge ids for.
  //!
  //! \param[in] e instance of entity to return wedge ids for.
  //!
  //! \return Return wedge ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) wedge_ids(E * e) const
  {
    return base_t::template entity_ids<wedge_t::dimension, wedge_t::domain>(e);
  }

  //============================================================================
  // Corner Interface
  //============================================================================

  //! \brief Return number of corners in the burton mesh.
  //! \return The number of corners in the burton mesh.
  size_t num_corners() const
  {
    return 
      base_t::template num_entities<corner_t::dimension, corner_t::domain>();
  }

  //! \brief Return all corners in the burton mesh.
  //!
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) corners() const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>();
  }

  //! \brief Return all corners in the burton mesh.
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) corners() // FIXME const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>();
  }


  //! \brief Return corners associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return corners for.
  //!
  //! \param[in] e instance of entity to return corners for.
  //!
  //! \return Return corners associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) corners(E * e) const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>(e);
  }

  //! \brief Return corners for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get corners for.
  //!
  //! \param[in] e Entity to get corners for.
  //!
  //! \return Corners for entity \e e in domain \e M.
  template<size_t M, class E>
  decltype(auto) corners(const flecsi::topology::domain_entity<M, E> & e) const
  {
    return 
      base_t::template entities<corner_t::dimension, M, corner_t::domain>(
        e.entity()
      );
  }

  //! \brief Return ids for all corners in the burton mesh.
  //! \return Ids for all corners in the burton mesh.
  decltype(auto) corner_ids() const
  {
    return base_t::template entity_ids<corner_t::dimension, corner_t::domain>();
  }

  //! \brief Return corner ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return corner ids for.
  //!
  //! \param[in] e instance of entity to return corner ids for.
  //!
  //! \return Return corner ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) corner_ids(E * e) const
  {
    return 
      base_t::template entity_ids<corner_t::dimension, corner_t::domain>(e);
  }


  //============================================================================
  // Region Interface
  //============================================================================

  //! \brief Return the number of regions in the burton mesh.
  //! \return The number of regions in the burton mesh.
  size_t num_regions() const
  {
    auto n = get_accessor(*this, mesh, num_regions, size_t, global, 0 );
    return *n;
  }

  //! \brief set the number of regions in the burton mesh.
  //! \param [in]  n  The number of regions in the burton mesh.
  void set_num_regions(size_t n)
  {
    auto n_acc = get_accessor(*this, mesh, num_regions, size_t, global, 0 ) ;
    *n_acc = n;
  }


  //! \brief Return the number of regions in the burton mesh.
  //! \param [in]  n  The number of regions in the burton mesh.
  template< typename T >
  void set_regions(T * region_ids)
  {
    for ( auto c : cells() )
      c->region() = region_ids[c.id()];
  }

  //! \brief Return all cells in the regions mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) regions()
  {
    // get the number of regions
    auto n = num_regions();
    // get the cells
    auto cs = cells();
    // create storage for regions
    using set_type_t = decltype( 
      base_t::template entities<cell_t::dimension, cell_t::domain>() 
    );

    // now filter out the cells of each region
    auto region_cells = cs.bin( [](const auto & c){ return c->region(); } );

    return region_cells;
  }

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
  }

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
  }

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
  }
  
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
    auto v = base_t::template make<vertex_t>( *this );
    base_t::template add_entity<vertex_t::dimension, vertex_t::domain>(v);
    v->coordinates() = pos;

    return v;
  }


  //============================================================================
  // Mesh Creation
  //============================================================================

  //! \brief Dump the burton mesh to standard out.
  std::ostream & dump( std::ostream & stream )
  {
    return base_t::dump( stream );
  }

  //! \brief Initialize burton mesh state for the number of \e vertices.
  //!
  //! \param[in] vertices The number of \e vertices to initialize the burton 
  //!                     mesh with.
  void init_parameters(size_t num_nodes)
  { }

  //!---------------------------------------------------------------------------
  //! \brief Initialize the burton mesh.
  //!---------------------------------------------------------------------------
  void init()
  {

    base_t::template init<0>();
    base_t::template init_bindings<1>();

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
        base_t::template reverse_entities<vertex_t::dimension, vertex_t::domain>(f);
      }
    } // for
#endif

    // register cell data
    register_data(*this, mesh, cell_volume, real_t, dense, 1, burton::cells);
    register_data(*this, mesh, cell_centroid, vector_t, dense, 1, burton::cells);
    register_data(*this, mesh, cell_min_length, real_t, dense, 1, burton::cells);

    // register face data
    register_data(*this, mesh, face_area, real_t, dense, 1, burton::faces);
    register_data(*this, mesh, face_normal, vector_t, dense, 1, burton::faces);
    register_data(*this, mesh, face_midpoint, vector_t, dense, 1, burton::faces);

    // register edge data
    register_data(*this, mesh, edge_midpoint, vector_t, dense, 1, burton::edges);
    
    // register wedge data
    register_data(*this, mesh, wedge_facet_area, real_t, dense, 1, burton::wedges);
    register_data(*this, mesh, wedge_facet_normal, vector_t, dense, 1, burton::wedges);
    register_data(*this, mesh, wedge_facet_centroid, vector_t, dense, 1, burton::wedges);
    
    // register time state
    register_data(*this, mesh, time, real_t, global, 1 );
    register_data(*this, mesh, time_step, size_t, global, 1 );

    auto soln_time = get_accessor(*this, mesh, time, real_t, global, 0);
    auto step = get_accessor(*this, mesh, time_step, size_t, global, 0);
    *soln_time = 0;
    *step = 0;

    // register some flags for identifying boundarys and various other things
    register_data(*this, mesh, node_flags, bitfield_t, dense, 1, burton::vertices);
    register_data(*this, mesh, edge_flags, bitfield_t, dense, 1, burton::edges);

    auto point_flags = get_accessor(*this, mesh, node_flags, bitfield_t, dense, 0);
    auto edge_flags = get_accessor(*this, mesh, edge_flags, bitfield_t, dense, 0);

    // register some flags for associating boundaries with entities
    register_data(*this, mesh, node_tags, tag_list_t, dense, 1, burton::vertices);
    register_data(*this, mesh, edge_tags, tag_list_t, dense, 1, burton::edges);
    register_data(*this, mesh, face_tags, tag_list_t, dense, 1, burton::faces);
    register_data(*this, mesh, cell_tags, tag_list_t, dense, 1, burton::cells);

    // now set the boundary flags.
    for ( auto f : faces() ) {
      // if there is only one cell, it is a boundary
      if ( f->is_boundary() ) {
        // point flags
        auto ps = vertices(f);
        for ( auto p : ps ) 
          point_flags[ p ].setbit( bits::boundary );
        // edge flags are only for 3d
        if ( num_dimensions == 3 ) {
          auto es = edges(f);
          for ( auto e : es ) 
            edge_flags[e].setbit( bits::boundary );
        } // dims
      } // is_boundary
    } // for

    // identify the cell regions
    register_data(*this, mesh, cell_region, size_t, dense, 1, burton::cells);
    register_data(*this, mesh, num_regions, size_t, global, 1);

    auto cell_region = get_accessor(*this, mesh, cell_region, size_t, dense, 0);
    auto num_regions = get_accessor(*this, mesh, num_regions, size_t, global, 0);

    *num_regions = 1;

    for ( auto c : cells() )
      cell_region[c] = 0;

    // update the geometry
    update_geometry();

  }


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


    //--------------------------------------------------------------------------
    // make sure face normal points out from first cell
    bool bad_face = false;
    auto fs = faces();
    auto num_faces = fs.size();

    #pragma omp parallel for reduction( || : bad_face )
    for( counter_t fid=0; fid<num_faces; ++fid ) {
      auto f = fs[fid];
      auto n = f->normal();
      auto fx = f->midpoint();
      auto c = cells(f).front();
      auto cx = c->midpoint();
      auto delta = fx - cx;
      auto dot = dot_product( n, delta );      
      if ( dot < 0 ) {
        bad_face = bad_face || true;
        #pragma omp critical
        {
          ss << "Face " << f.id() << " has opposite normal" << std::endl;
        }
      }
    } 

    if ( bad_face ) return raise_or_return( ss );


    //--------------------------------------------------------------------------
    // check all the corners and wedges
    bool bad_corner = false;
    auto cnrs = corners();
    auto num_corners = cnrs.size();

    //#omp parallel for reduction( || : bad_corner )
    for( counter_t cnid=0; cnid<num_corners; ++cnid ) {

      auto cn = cnrs[cnid];
      auto cs = cells(cn);
      auto fs = faces(cn);
      auto es = edges(cn);
      auto vs = vertices(cn);
      auto ws = wedges(cn);

      if ( cs.size() != 1 ) {
        #pragma omp critical
        ss << "Corner " << cn.id() << " has " << cs.size() << "/=1 cells" 
           << std::endl;
      }

      if ( fs.size() != num_dimensions ) {
        #pragma omp critical
        {
          ss << "Corner " << cn.id() << " has " << fs.size() << "/=" 
             << num_dimensions << " faces" << std::endl;
        }
      }

      if ( es.size() != num_dimensions ) {
        #pragma omp critical
        {
          ss << "Corner " << cn.id() << " has " << es.size() << "/=" 
             << num_dimensions << " edges" << std::endl;
        }
      }

      if ( vs.size() != 1 ) {
        #pragma omp critical
        {
          ss << "Corner " << cn.id() << " has " << vs.size() << "/=1 vertices"
             << std::endl;
        }
      }

      auto cl = cs.front();
      auto vt = vs.front();
      
      if ( ws.size() % 2 != 0 ) {
        #pragma omp critical
        {
          ss << "Corner " << cn.id() << " has " << ws.size() << "%2/=0 wedges"
             << std::endl;
        }
      }

      for ( auto wg = ws.begin(); wg != ws.end();  ) 
        for ( auto i=0; i<2 && wg != ws.end(); i++, ++wg)
        {
          auto cls = cells( *wg );
          auto fs = faces( *wg );
          auto es = edges( *wg );
          auto vs = vertices( *wg );
          auto cns = corners( *wg );
          if ( cls.size() != 1 ) {
            #pragma omp critical
            {          
              ss << "Wedge " << (*wg).id() << " has " << cls.size() 
                 << "/=1 cells" << std::endl;
            }
          }
          if ( fs.size() != 1 ) {
            #pragma omp critical 
            {           
              ss << "Wedge " << (*wg).id() << " has " << fs.size() 
                 << "/=1 faces" << std::endl;
            }
          }
          if ( es.size() != 1 ) {
            #pragma omp critical 
            {
              ss << "Wedge " << (*wg).id() << " has " << es.size() 
                 << "/=1 edges" << std::endl;
            }
          }
          if ( vs.size() != 1 ) {
            #pragma omp critical
            {
              ss << "Wedge " << (*wg).id() << " has " << vs.size() 
                 << "/=1 vertices" << std::endl;
            }
          }
          if ( cns.size() != 1 ) {
            #pragma omp critical
            {
              ss << "Wedge " << (*wg).id() << " has " << cns.size() 
                 << "/=1 corners" << std::endl;
            }
          }
          auto vert = vs.front();
          auto cell = cls.front();
          auto corn = cns.front();
          if ( vert != vt ) {
            #pragma omp critical
            {
              ss << "Wedge " << (*wg).id() << " has incorrect vertex " 
                 << vert.id() << "!=" << vt.id() << std::endl;
            }
          }
          if ( cell != cl ) {
            #pragma omp critical
            {
              ss << "Wedge " << (*wg).id() << " has incorrect cell " 
                 << cell.id() << "!=" << cl.id() << std::endl;
            }
          }
          if ( corn != cn ) {
            #pragma omp critical
            {
              ss << "Wedge " << (*wg).id() << " has incorrect corner " 
                 << corn.id() << "!=" << cn.id() << std::endl;
            }
          }
          auto fc = fs.front();            
          auto fx = fc->midpoint();
          auto cx = cl->midpoint();
          auto delta = fx - cx;
          real_t dot;
          if ( i == 0 ) {
            auto n = (*wg)->facet_normal_right();
            dot = dot_product( n, delta );
          }
          else {
            auto n = (*wg)->facet_normal_left();
            dot = dot_product( n, delta );
          }
          if ( dot < 0 ) {
            #pragma omp critical
            {
              ss << "Wedge " << (*wg).id() << " has opposite normal" 
                 << std::endl;
            }
          }
        } // wedges
      
    } // corners


    if ( bad_corner ) return raise_or_return( ss );


    return true;

  }

  //!---------------------------------------------------------------------------
  //! \brief Compute the goemetry.
  //!---------------------------------------------------------------------------
  void update_geometry()
  {
    // get the mesh info
    auto cs = cells();
    auto fs = faces();
    auto es = edges();  
    auto cnrs = corners();
    auto num_cells = cs.size();
    auto num_faces = fs.size();
    auto num_edges = es.size();
    auto num_corners = cnrs.size();

    // get all the data now so we can put everything in one parallel region
    auto cell_center = get_accessor(*this, mesh, cell_centroid, vector_t, dense, 0);
    auto cell_volume = get_accessor(*this, mesh, cell_volume, real_t, dense, 0);
    auto cell_min_length = get_accessor(*this, mesh, cell_min_length, real_t, dense, 0);

    auto face_area = get_accessor(*this, mesh, face_area, real_t, dense, 0);
    auto face_norm = get_accessor(*this, mesh, face_normal, vector_t, dense, 0);
    auto face_midp = get_accessor(*this, mesh, face_midpoint, vector_t, dense, 0); 

    auto edge_midp = get_accessor(*this, mesh, edge_midpoint, vector_t, dense, 0); 

    auto wedge_facet_normal = get_accessor(*this, mesh, wedge_facet_normal, vector_t, dense, 0);
    auto wedge_facet_area = get_accessor(*this, mesh, wedge_facet_area, real_t, dense, 0);
    auto wedge_facet_centroid = get_accessor(*this, mesh, wedge_facet_centroid, vector_t, dense, 0); 

    //--------------------------------------------------------------------------
    // compute cell parameters

    #pragma omp parallel
    {

      #pragma omp for   nowait
      for ( counter_t i=0; i<num_cells; i++ ) {
        auto c = cs[i];
        cell_volume[c] = c->volume();
        cell_center[c] = c->centroid();
        cell_min_length[c] = c->min_length();
      } 

      //--------------------------------------------------------------------------
      // compute face parameters

      #pragma omp for nowait
      for ( counter_t i=0; i<num_faces; i++ ) {
        auto f = fs[i];
        face_area[f] = f->area();
        face_norm[f] = f->normal() / face_area[f];
        face_midp[f] = f->midpoint();
      } 

      //--------------------------------------------------------------------------
      // compute edge parameters

      #pragma omp for
      for ( counter_t i=0; i<num_edges; i++ ) {
        auto e = es[i];
        edge_midp[e] = e->midpoint();
      } 

      //--------------------------------------------------------------------------
      // compute wedge parameters

      #pragma omp for
      for ( counter_t i=0; i<num_corners; ++i ) {
        auto cn = cnrs[i];
        auto ws = wedges(cn);
        // both wedges SHOULD have the same vertex
        const auto & v = vertices(cn).front()->coordinates();
        // first compute the normals
        for ( auto wit = ws.begin(); wit != ws.end(); ++wit ) 
        {
          // get the first wedge normal
          {
            const auto & e = edge_midp[ edges(*wit).front() ];
            const auto & f = face_midp[ faces(*wit).front() ];
            wedge_facet_normal[*wit] = wedge_t::facet_normal_right( v, e, f );
          }
          // move to next wedge
          ++wit;
          assert( wit != ws.end() );
          // get the second wedge normal
          {
            const auto & e = edge_midp[ edges(*wit).front() ];
            const auto & f = face_midp[ faces(*wit).front() ];
            wedge_facet_normal[*wit] = wedge_t::facet_normal_left( v, e, f );
          }
        }
        // now normalize the normals and compute other quantities
        for ( auto w : ws) {
          wedge_facet_area[w] = abs( wedge_facet_normal[w] );
          wedge_facet_normal[w] /= wedge_facet_area[w];
          wedge_facet_centroid[w] = w->facet_centroid();
        }
      }

    } // end omp parallel

  }


  //============================================================================
  //! \brief Install a boundary and tag the relatex entities.
  //============================================================================
  template< typename P >
  tag_t install_boundary( P && p ) 
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
        // tag the face
        f->tag( this_bnd );
        this_bnd_faces.emplace_back( f );
        // tag the vertices
        auto vs = vertices( f );
        this_bnd_verts.reserve( this_bnd_verts.size() + vs.size() );
        for ( auto v : vs ) this_bnd_verts.emplace_back( v );
        // tag edges in 3d
        if ( num_dimensions == 3 ) {
          auto es = edges( f );
          this_bnd_edges.reserve( this_bnd_edges.size() + es.size() );
          for ( auto e : es ) this_bnd_edges.emplace_back( e );
        } // dims
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
    for ( auto e : this_bnd_edges ) e->tag( this_bnd );
    // add the vertex tags
    for ( auto v : this_bnd_verts ) v->tag( this_bnd );

    // if it's two dimensions, copy the face list into the edge list
    // if ( num_dimensions == 2 )
    //   this_bnd_edges.assign( this_bnd_faces.begin(), this_bnd_faces.end() );

    return this_bnd;
  }

  //============================================================================
  //! \brief Get the set of tagged vertices associated with a specific id
  //! \praram [in] id  The tag to lookup.
  //! \return The set of tagged vertices.
  //============================================================================
  const auto & tagged_vertices( tag_t id ) const noexcept
  {
    return vert_sets_[ id ];
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
    
    E * e = nullptr;

    switch ( verts.size() ) {
    case (1,2):
      raise_runtime_error( "can't have <3 vertices" );
    case (3):
      e = base_t::template make< burton_triangle_t<num_dimensions> >(*this);
      break;
    case (4):
      e = base_t::template make< burton_quadrilateral_t<num_dimensions> >(*this);
      break;
    default:
      e = base_t::template make< burton_polygon_t<num_dimensions> >(*this);
      break;
    }

    base_t::template add_entity<E::dimension, E::domain>( e );
    base_t::template init_entity<E::domain, E::dimension, vertex_t::dimension>( e, std::forward<V>(verts) );
    return e;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_3d_element_from_verts_( V && verts )
  {
    
    cell_t * c = nullptr;

    switch ( verts.size() ) {
    case (1,2,3):
      raise_runtime_error( "can't have <4 vertices" );
    case (4):
      c = base_t::template make< burton_tetrahedron_t >(*this);
      break;
    case (8):
      c = base_t::template make< burton_hexahedron_t >(*this);
      break;
    default:
      raise_runtime_error( "can't build polyhedron from vertices alone" );      
      break;
    }

    base_t::template add_entity<cell_t::dimension, cell_t::domain>(c);
    base_t::template init_cell<cell_t::domain>( c, std::forward<V>(verts) );
    return c;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  template< typename F >
  auto create_3d_element_from_faces_( F && faces )
  {
    
    cell_t * c = nullptr;

    switch ( faces.size() ) {
    case (1,2,3):
      raise_runtime_error( "can't have <4 vertices" );
    default:
      c = base_t::template make< burton_polyhedron_t >( *this, std::forward<F>(faces) );
      break;
    }

    base_t::template add_entity<cell_t::dimension, cell_t::domain>(c);
    base_t::template init_entity<cell_t::domain, cell_t::dimension, face_t::dimension>( c, std::forward<F>(faces) );
    return c;
  } // create_cell


  //============================================================================
  // Private Data 
  //============================================================================

  //! \brief Tagged sets 
  //@ {
  std::vector< std::vector<face_t*> >   face_sets_;
  std::vector< std::vector<edge_t*> >   edge_sets_;
  std::vector< std::vector<vertex_t*> > vert_sets_;
  //@ }


}; // class burton_mesh_t


////////////////////////////////////////////////////////////////////////////////
// External Class Definitions
////////////////////////////////////////////////////////////////////////////////


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
