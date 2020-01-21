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
#include <flecsale/eqns/euler_eqns.h>
#include <flecsale/eqns/flux.h>
#include <flecsale/eos/ideal_gas.h>
#include <ristra/math/general.h>

#include <flecsi-sp/utils/char_array.h>
#include <flecsi-sp/utils/types.h>
#include <flecsi-sp/burton/burton_mesh.h>

#include <flecsi/data/global_accessor.h>

#include "../common/utils.h"

namespace apps {
namespace hydro {

// mesh and some underlying data types
using mesh_t = flecsi_sp::burton::burton_mesh_t;
using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;
using counter_t = mesh_t::counter_t;

using eos_t = flecsale::eos::ideal_gas_t<real_t>;

using eqns_t = typename flecsale::eqns::euler_eqns_t<real_t, mesh_t::num_dimensions>;

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

template<typename T>
using global_handle_w = flecsi::global_accessor_u<T, flecsi::wo>;

template<typename T>
using global_handle_rw = flecsi::global_accessor_u<T, flecsi::rw>;

template<typename T>
using global_handle_r = flecsi::global_accessor_u<T, flecsi::ro>;

template<typename DC>
using client_handle_w = flecsi_sp::utils::client_handle_w<DC>;

template<typename DC>
using client_handle_r = flecsi_sp::utils::client_handle_r<DC>;

template<typename T>
using future_handle =
  flecsi::execution::flecsi_future<
    T, flecsi::execution::launch_type_t::index
  >;


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
//! \brief alias the flux function
//! Change the called function to alter the flux evaluation.
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename UL, typename UR, typename V >
FLECSI_INLINE_TARGET auto flux_function( UL && left_state, UR && right_state, V && norm )
{ 
  return 
    flecsale::eqns::hlle_flux<E>(
        std::forward<UL>(left_state), 
        std::forward<UR>(right_state), 
        std::forward<V>(norm) ); 
}

////////////////////////////////////////////////////////////////////////////////
//! \brief alias the boundary flux function
//! Change the called function to alter the flux evaluation.
////////////////////////////////////////////////////////////////////////////////
template< typename E, typename U, typename V >
FLECSI_INLINE_TARGET auto boundary_flux( U && state, V && norm )
{ 
  return 
    E::wall_flux( std::forward<U>(state), std::forward<V>(norm) ); 
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Pack data into a tuple
//! Change the called function to alter the flux evaluation.
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename...ARGS >
FLECSI_INLINE_TARGET decltype(auto) pack( T && loc, ARGS&&... args )
{ 
  return 
    std::forward_as_tuple( std::forward<ARGS>(args)(std::forward<T>(loc))... ); 
}

} // namespace hydro
} // namespace apps
