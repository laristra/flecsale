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
#include <flecsale-config.h>
#include <flecsale/eqns/lagrange_eqns.h>
#include <flecsale/eqns/flux.h>
#include <flecsale/eos/ideal_gas.h>
#include <ristra/math/general.h>
#include <ristra/math/matrix.h>

#include <flecsi-sp/utils/char_array.h>
#include <flecsi-sp/utils/types.h>
#include <flecsi-sp/burton/burton_mesh.h>

#include "../common/utils.h"

namespace apps {
namespace hydro {

// mesh and some underlying data types
using mesh_t = flecsi_sp::burton::burton_mesh_t;
using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;
using counter_t = mesh_t::counter_t;

using eos_t = flecsale::eos::ideal_gas_t<real_t>;

using eqns_t = typename flecsale::eqns::lagrange_eqns_t<real_t, mesh_t::num_dimensions>;

template< std::size_t N >
using matrix = ristra::math::matrix< real_t, N, N >; 

using flux_data_t = eqns_t::flux_data_t;


// explicitly use some other stuff
using std::cout;
using std::cerr;
using std::endl;


// the access permission types
template<typename T>
using dense_handle_w = flecsi_sp::utils::dense_handle_w<T>;

template<typename T>
using dense_handle_rw = flecsi_sp::utils::dense_handle_rw<T>;

template<typename T>
using dense_handle_r = flecsi_sp::utils::dense_handle_r<T>;

template<typename DC>
using client_handle_w = flecsi_sp::utils::client_handle_w<DC>;

template<typename DC>
using client_handle_r = flecsi_sp::utils::client_handle_r<DC>;

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
using char_array_t = flecsi_sp::utils::char_array_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief A general boundary condition type.
//! \tparam N  The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
class boundary_condition_t
{
public:

  using real_type   = real_t; 
  using vector_type = vector_t; 

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
class symmetry_boundary_condition_t : public boundary_condition_t
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
inline boundary_condition_t * make_boundary_condition( const std::string & str )
{
  if ( str == "symmetry" )
    return new symmetry_boundary_condition_t();
  else if ( str == "none" )
    return new boundary_condition_t();
  else {
    throw_implemented_error( 
      "No implementation for boundary condition of type \'" << str << "\'"
    );
    return nullptr;
  }

}


//! \brief a type for storing boundary tags
using tag_t = mesh_t::tag_t;

//! \brief a map for storing links between boundary conditions and tags
using boundary_map_t = std::map< tag_t, boundary_condition_t * >;

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
