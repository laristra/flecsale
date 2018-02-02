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

#include "../common/utils.h"

namespace apps {
namespace hydro {

// mesh and some underlying data types
using mesh_t = flecsi_sp::burton::burton_mesh_t;
using real_t = mesh_t::real_t;
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
using dense_handle_w__ = flecsi_sp::utils::dense_handle_w__<T>;

template<typename T>
using dense_handle_rw__ = flecsi_sp::utils::dense_handle_rw__<T>;

template<typename T>
using dense_handle_r__ = flecsi_sp::utils::dense_handle_r__<T>;

template<typename DC>
using client_handle_w__ = flecsi_sp::utils::client_handle_w__<DC>;

template<typename DC>
using client_handle_r__ = flecsi_sp::utils::client_handle_r__<DC>;

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
auto flux_function( UL && left_state, UR && right_state, V && norm )
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
auto boundary_flux( U && state, V && norm )
{ 
  return 
    E::wall_flux( std::forward<U>(state), std::forward<V>(norm) ); 
}

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

} // namespace hydro
} // namespace apps
