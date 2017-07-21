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
#include <flecsale/eqns/euler_eqns.h>
#include <flecsale/eqns/flux.h>
#include <flecsale/eos/ideal_gas.h>
#include <flecsale/math/general.h>

#include <flecsale/mesh/burton/burton.h>


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
//namespace io    = flecsale::io;

// the handle type
template<typename T, size_t EP, size_t SP, size_t GP>
using dense_handle__ =
  flecsi::data::legion::dense_handle_t<T, EP, SP, GP>;

template<typename DC, size_t PS>
using client_handle__ = flecsi::data_client_handle__<DC, PS>;

// mesh and some underlying data types
template <std::size_t N>
using mesh__ = typename mesh::burton::burton_mesh_t<N>;

using size_t = common::size_t;
using real_t = common::real_t;

using eos_t = eos::eos_base_t<real_t>;

template< std::size_t N >
using eqns__ = typename eqns::euler_eqns_t<real_t, N>;

template< std::size_t N >
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

////////////////////////////////////////////////////////////////////////////////
//! \brief alias the flux function
//! Change the called function to alter the flux evaluation.
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename UL, typename UR, typename V >
auto flux_function( UL && left_state, UR && right_state, V && norm )
{ 
  return 
    eqns::hlle_flux<E>( std::forward<UL>(left_state), 
                        std::forward<UR>(right_state), 
                        std::forward<V>(norm) ); 
}

////////////////////////////////////////////////////////////////////////////////
//! \brief alias the boundary flux function
//! Change the called function to alter the flux evaluation.
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
auto boundary_flux( U && state, V && norm )
{ 
  return 
    E::wall_flux( std::forward<U>(state), std::forward<V>(norm) ); 
}

} // namespace hydro
} // namespace apps
