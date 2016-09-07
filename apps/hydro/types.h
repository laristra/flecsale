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
#include <ale/common/types.h>
#include <ale/eqns/euler_eqns.h>
#include <ale/eqns/flux.h>
#include <ale/eos/ideal_gas.h>
#include <ale/math/general.h>

#include <ale/mesh/burton/burton.h>


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

using eos_t = eos::ideal_gas_t<real_t>;

template< std::size_t N >
using eqns_t = typename eqns::euler_eqns_t<real_t, N>;

template< std::size_t N >
using flux_data_t = typename eqns_t<N>::flux_data_t;

// explicitly use some other stuff
using std::cout;
using std::cerr;
using std::endl;

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


////////////////////////////////////////////////////////////////////////////////
//! \brief A functor for accessing state in the mesh
//! \tparam M  the mesh type
////////////////////////////////////////////////////////////////////////////////
template < typename M >
class state_accessor 
{
public:

  //! typedefs
  using real_t = typename M::real_t;
  using vector_t = typename M::vector_t;

  //! \brief determine the type of accessor
  //! \tparam T the data type we are accessing
  template< typename T >
  using accessor_t = 
    decltype( access_state( std::declval<M>(), "var_name", T ) );

  //! \brief main constructor
  //! \param [in]  mesh  the mesh object
  constexpr state_accessor( M & mesh ) :
    d( access_state( mesh, "density",           real_t ) ),
    p( access_state( mesh, "pressure",          real_t ) ),
    v( access_state( mesh, "velocity",        vector_t ) ),
    e( access_state( mesh, "internal_energy",   real_t ) ),
    t( access_state( mesh, "temperature",       real_t ) ),
    a( access_state( mesh, "sound_speed",       real_t ) )
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
    return std::forward_as_tuple( d[ forward<T>(i) ], 
                                  v[ forward<T>(i) ], 
                                  p[ forward<T>(i) ], 
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
    return std::forward_as_tuple( d[ forward<T>(i) ], 
                                  v[ forward<T>(i) ], 
                                  p[ forward<T>(i) ], 
                                  e[ forward<T>(i) ], 
                                  t[ forward<T>(i) ], 
                                  a[ forward<T>(i) ] );
  }

private:

  // private data
  accessor_t<real_t>   d;
  accessor_t<real_t>   p;
  accessor_t<vector_t> v;
  accessor_t<real_t>   e;
  accessor_t<real_t>   t;
  accessor_t<real_t>   a;
       
};

} // namespace hydro
} // namespace apps
