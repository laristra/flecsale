/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Define the main types for the hydro solver.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsale/common/types.h>
#include <flecsale/eqns/lagrange_eqns.h>
#include <flecsale/eqns/flux.h>
#include <flecsale/eos/ideal_gas.h>
#include <flecsale/math/general.h>
#include <flecsale/math/matrix.h>

#include <flecsale/mesh/burton/burton.h>


namespace apps {
namespace hydro {

// some namespace aliases
namespace common= ale::common;
namespace mesh  = ale::mesh;
namespace math  = ale::math;
namespace utils = ale::utils;
namespace geom  = ale::geom;
namespace eos   = ale::eos;
namespace eqns  = ale::eqns;

// mesh and some underlying data types
using mesh_2d_t = mesh::burton_mesh_2d_t;
using mesh_3d_t = mesh::burton_mesh_3d_t;

using size_t = common::size_t;
using real_t = common::real_t;

template< std::size_t N >
using matrix_t = math::matrix< real_t, N, N >; 

using eos_t = eos::eos_base_t<real_t>;

template<std::size_t N>
using eqns_t = eqns::lagrange_eqns_t<real_t, N>;

template<std::size_t N>
using flux_data_t = typename eqns_t<N>::flux_data_t;

// explicitly use some other stuff
using std::cout;
using std::cerr;
using std::endl;

//! \brief a class to distinguish between different types of 
//!   update errors.
enum class solution_error_t 
{
  unphysical, variance, ok
};

//! \brief a class to distinguish between different types of 
//!   restarting modes.
enum class mode_t 
{
  normal, retry, restart, quit
};


////////////////////////////////////////////////////////////////////////////////
//! \brief A general boundary condition type.
//! \tparam N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class boundary_condition_t
{
public:

  using real_type   = real_t; 
  using vector_type = math::vector< real_type, N >; 

  virtual bool has_prescribed_velocity() const 
  { return false; };

  virtual bool has_prescribed_pressure() const 
  { return false; };

  virtual bool has_symmetry() const 
  { return false; };

  virtual vector_type velocity( const vector_type & x, const real_type & t ) const 
  { return 0; };

  virtual real_type pressure( const vector_type & x, const real_type & t ) const 
  { return 0; };

  virtual ~boundary_condition_t() {}
};

////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the boundary condition type for symmetry 
//!        conditions.
//! \tparam N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class symmetry_boundary_condition_t : public boundary_condition_t<N>
{
public:
  virtual bool has_symmetry() const override
  { return true; }

  virtual ~symmetry_boundary_condition_t() {}
};

////////////////////////////////////////////////////////////////////////////////
//! \brief A function to create a new boundary object based on a string.
//! \tparam N  The number of dimensions.
//! \param [in] type_str  The type of boundary condition as a string.
//! \return A new instance of the type.
////////////////////////////////////////////////////////////////////////////////
template < std::size_t N >
boundary_condition_t<N> * make_boundary_condition( const std::string & str )
{
  if ( str == "symmetry" )
    return new symmetry_boundary_condition_t<N>();
  else if ( str == "none" )
    return new boundary_condition_t<N>();
  else {
    raise_implemented_error( 
      "No implementation for boundary condition of type \'" << str << "\'"
    );
    return nullptr;
  }

}


//! \brief a type for storing boundary tags
using tag_t = mesh_2d_t::tag_t;

//! \brief a map for storing links between boundary conditions and tags
template< std::size_t N >
using boundary_map_t = std::map< tag_t, boundary_condition_t<N> * >;

//! \breif a map for equations of state
using eos_map_t = std::map< tag_t, eos_t * >;

////////////////////////////////////////////////////////////////////////////////
//! \brief A functor for accessing state in the mesh
//! \tparam M  the mesh type
////////////////////////////////////////////////////////////////////////////////
template < typename M >
class cell_state_accessor 
{
public:
  
  //! some type aliases
  using real_t = typename M::real_t;
  using vector_t = typename M::vector_t;

  //! \brief determine the type of accessor
  //! \tparam T the data type we are accessing
  template< typename T >
  using accessor_t = 
    decltype( flecsi_get_accessor( std::declval<M>(), hydro, var_name, T, dense, 0 ) );

  //! \brief main constructor
  //! \param [in]  mesh  the mesh object
  constexpr cell_state_accessor( M & mesh ) :
    m( flecsi_get_accessor( mesh, hydro, cell_mass, real_t, dense, 0 ) ),
    V( flecsi_get_accessor( mesh, hydro, cell_volume, real_t, dense, 0 ) ),
    p( flecsi_get_accessor( mesh, hydro, cell_pressure, real_t, dense, 0 ) ),
    v( flecsi_get_accessor( mesh, hydro, cell_velocity, vector_t, dense, 0 ) ),
    d( flecsi_get_accessor( mesh, hydro, cell_density, real_t, dense, 0 ) ),
    e( flecsi_get_accessor( mesh, hydro, cell_internal_energy, real_t, dense, 0 ) ),
    t( flecsi_get_accessor( mesh, hydro, cell_temperature, real_t, dense, 0 ) ),
    a( flecsi_get_accessor( mesh, hydro, cell_sound_speed, real_t, dense, 0 ) )
  {}

  //! \brief main accessor
  //! \remark const version
  //!
  //! \tparam T the element type
  //! \param [in] i  the element we are accessing
  //! \return a tuple of references to all the state data
  template< typename T >
  constexpr auto operator()( T && i ) const 
  {
    using std::forward;
    return std::forward_as_tuple( V[ forward<T>(i) ], 
                                  m[ forward<T>(i) ], 
                                  v[ forward<T>(i) ], 
                                  p[ forward<T>(i) ], 
                                  d[ forward<T>(i) ], 
                                  e[ forward<T>(i) ], 
                                  t[ forward<T>(i) ], 
                                  a[ forward<T>(i) ] );
  }

  //! \brief main accessor
  //! \remark non-const version
  //!
  //! \tparam T the element type
  //! \param [in] i  the element we are accessing
  //! \return a tuple of references to all the state data
  template< typename T >
  auto operator()( T && i ) 
  {
    using std::forward;
    return std::forward_as_tuple( V[ forward<T>(i) ], 
                                  m[ forward<T>(i) ], 
                                  v[ forward<T>(i) ], 
                                  p[ forward<T>(i) ], 
                                  d[ forward<T>(i) ], 
                                  e[ forward<T>(i) ], 
                                  t[ forward<T>(i) ], 
                                  a[ forward<T>(i) ] );
  }

private:

  // private data
  accessor_t<real_t>   m;
  accessor_t<real_t>   V;
  accessor_t<real_t>   p;
  accessor_t<vector_t> v;
  accessor_t<real_t>   d;
  accessor_t<real_t>   e;
  accessor_t<real_t>   t;
  accessor_t<real_t>   a;
       
};

////////////////////////////////////////////////////////////////////////////////
//! \brief A class to make setting the cfl easy
////////////////////////////////////////////////////////////////////////////////
struct time_constants_t {

  real_t accoustic = 1.0;
  real_t volume = 1.0;
  real_t growth = 0.0; 

};

} // namespace hydro
} // namespace apps
