/*~-------------------------------------------------------------------------~~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
/*!
 *
 * \file types.h
 * 
 * \brief Define the main types for the hydro solver.
 *
 ******************************************************************************/

#pragma once

// user includes
#include <ale/eqns/lagrange_eqns.h>
#include <ale/eqns/flux.h>
#include <ale/eos/ideal_gas.h>
#include <ale/math/math.h>
#include <ale/math/matrix.h>

#include <ale/mesh/burton/burton.h>


namespace apps {
namespace hydro {

// some namespace aliases
namespace mesh  = ale::mesh;
namespace math  = ale::math;
namespace utils = ale::utils;
namespace geom  = ale::geom;
namespace eos   = ale::eos;
namespace eqns  = ale::eqns;

// mesh and some underlying data types
using mesh_t = mesh::burton_mesh_t;

using size_t = mesh_t::size_t;
using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;
using matrix_t = math::matrix< real_t, mesh_t::num_dimensions(), mesh_t::num_dimensions() >;

using eos_t = eos::ideal_gas_t<real_t>;

using eqns_t = eqns::lagrange_eqns_t<real_t, mesh_t::num_dimensions()>;
using flux_data_t = eqns_t::flux_data_t;

// explicitly use some other stuff
using std::cout;
using std::cerr;
using std::endl;



////////////////////////////////////////////////////////////////////////////////
//! \brief A functor for accessing state in the mesh
//! \tparam M  the mesh type
////////////////////////////////////////////////////////////////////////////////
template < typename M >
class cell_state_accessor 
{
public:

  //! \brief determine the type of accessor
  //! \tparam T the data type we are accessing
  template< typename T >
  using accessor_t = 
    decltype( access_state( std::declval<M>(), "var_name", T ) );

  //! \brief main constructor
  //! \param [in]  mesh  the mesh object
  constexpr cell_state_accessor( M & mesh ) :
    m( access_state( mesh, "cell_mass",              real_t ) ),
    V( access_state( mesh, "cell_volume",            real_t ) ),
    p( access_state( mesh, "cell_pressure",          real_t ) ),
    v( access_state( mesh, "cell_velocity",        vector_t ) ),
    d( access_state( mesh, "cell_density",           real_t ) ),
    e( access_state( mesh, "cell_internal_energy",   real_t ) ),
    t( access_state( mesh, "cell_temperature",       real_t ) ),
    a( access_state( mesh, "cell_sound_speed",       real_t ) )
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
