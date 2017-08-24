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
#include <flecsale/utils/trivial_string.h>

#include <flecsale/mesh/burton/burton.h>

#include "../common/utils.h"

namespace apps {
namespace hydro {

// some namespace aliases
namespace common= flecsale::common;
namespace mesh  = flecsale::mesh;
namespace math  = flecsale::math;
namespace utils = flecsale::utils;
namespace geom  = flecsale::geom;
namespace eos   = flecsale::eos;
namespace eqns  = flecsale::eqns;

// the handle type
template<typename T>
using dense_handle_w__ =
  flecsi::data::legion::dense_handle_t<T, flecsi::wo, flecsi::wo, flecsi::ro>;

template<typename T>
using dense_handle_rw__ =
  flecsi::data::legion::dense_handle_t<T, flecsi::rw, flecsi::rw, flecsi::ro>;

template<typename T>
using dense_handle_r__ =
  flecsi::data::legion::dense_handle_t<T, flecsi::ro, flecsi::ro, flecsi::ro>;

template<typename DC>
using client_handle_w__ = flecsi::data_client_handle__<DC, flecsi::wo>;

template<typename DC>
using client_handle_r__ = flecsi::data_client_handle__<DC, flecsi::ro>;


// mesh and some underlying data types
template <std::size_t N>
using mesh__ = typename mesh::burton::burton_mesh_t<N, true>;

using size_t = common::size_t;
using real_t = common::real_t;

template< std::size_t N >
using matrix__ = math::matrix< real_t, N, N >; 

using eos_t = eos::ideal_gas_t<real_t>;

template< std::size_t N >
using eqns__ = typename eqns::lagrange_eqns_t<real_t, N>;

template<std::size_t N>
using flux_data__ = typename eqns__<N>::flux_data_t;

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

//! a trivially copyable character array
using char_array_t = utils::trivial_string_t;

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
using tag_t = mesh__<2>::tag_t;

//! \brief a map for storing links between boundary conditions and tags
template< std::size_t N >
using boundary_map_t = std::map< tag_t, boundary_condition_t<N> * >;

//! \breif a map for equations of state
using eos_map_t = std::map< tag_t, eos_t * >;

////////////////////////////////////////////////////////////////////////////////
//! \brief Pack data into a tuple
//! Change the called function to alter the flux evaluation.
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename...ARGS >
decltype(auto) pack( T && loc, ARGS&&... args )
{ 
  return 
    std::forward_as_tuple( std::forward<ARGS>(args)(std::forward<T>(loc))... ); 
}


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
