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
 * \file lagrange_eqns.h
 * 
 * \brief The desrciption of the euler equations.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/math/tuple.h"
#include <ale/math/math.h>
#include "ale/math/vector.h"

namespace ale {
namespace eqns {

////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the euler equations
////////////////////////////////////////////////////////////////////////////////
template<typename T, size_t N>
struct lagrange_eqns_t {


public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the size type
  using size_t = std::size_t;

  // \brief the real type
  using real_t = T;

  //! \brief the vector_type
  using vector_t = math::vector<real_t,N>;

  //! the number of dimensions
  static constexpr size_t dimensions = N;

  //============================================================================
  // \brief The equations struct
  //============================================================================
  struct equations {

    //! \brief  the type for holding the state data (mass, momentum, and energy)
    //! 
    //! tuple_t is a std::tuple with a real_t for mass, a vector_t for momentum
    //! and a real_t for energy.  This needs to correspond to index
    //! or there may be problems
    using data_t = math::tuple<real_t,vector_t,real_t>;

    //! \brief the variables in the primitive state
    enum index : size_t { 
      mass = 0, 
      momentum, 
      energy,
      total
    };

    //! \brief the number of equations
    static constexpr size_t number(void)
    {  return index::total; }

  };

  //============================================================================
  // \brief The variables struct
  //============================================================================
  struct variables {


    //! \brief  the type for holding the state data
    //!
    //! This needs to correspond with index or they may be problems
    using data_t = 
      math::tuple<real_t,real_t,vector_t,real_t,real_t,real_t,real_t,real_t>;

    //! \brief the variables in the primitive state
    enum index : size_t
    { 
      volume = 0, 
      mass,
      velocity, 
      pressure,
      density, 
      internal_energy,
      temperature,
      sound_speed,
      total
    };

    //! \brief the number of variables
    static constexpr size_t number(void)
    {  return index::total; }
    
  };


  //============================================================================
  // Deferred Typedefs
  //============================================================================

  //! \brief  the type for holding the state data
  using state_data_t = typename variables::data_t;
  //! \brief  the type for holding refernces to state data
  using state_ref_t = utils::reference_wrapper_t< state_data_t >;

  //! \brief  the type for holding the state data (mass, momentum, and energy)
  using flux_data_t = typename equations::data_t;
  //! \brief  the type for holding refernces to flux data
  using flux_ref_t = utils::reference_wrapper_t< flux_data_t >;



  //============================================================================
  //! \brief accessors for various quantities
  //!
  //! Wrapper functions are only used for accessing either independant or 
  //! derived quantities.  Modifictions of variables should access
  //! the element of the state_data_t explicitily, because it should remain 
  //! clear exactly what variable is being set.  Does 
  //!   set_density(val) 
  //! set density explicitly? Or does it really set temperature and pressure?
  //! But it is clear that
  //!   get<variables::index::density> = val
  //! is assigning a new value to density.
  //!
  //! \param [in] u the state
  //! \return the quantity of interest
  //============================================================================
  static auto volume( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::volume>( u ); 
  }

  static auto mass( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::mass>( u ); 
  }

  static auto velocity( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::velocity>( u ); 
  }

  static auto pressure(  const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::pressure>( u ); 
  }

  static auto density( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::density>( u ); 
  }

  static auto internal_energy( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::internal_energy>( u ); 
  }

  static auto sound_speed( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::sound_speed>( u ); 
  }

  static auto total_energy( const state_data_t & u )
  { 
    using math::abs;
    auto ie = internal_energy( u );
    auto vel = velocity( u );
    return ie + 0.5 * abs( vel );
  }

  static auto impedance( const state_data_t & u )
  { 
    return density(u) * sound_speed(u);
  }

  //============================================================================
  //! \brief Compute the fastest moving eigenvalue
  //! \param [in]  u   The solution state
  //! \return the fastest moving wave speed
  //============================================================================
  static auto fastest_wavespeed( const state_data_t & u )
  {
    auto a = sound_speed( u );
    return a;
  }

  //============================================================================
  //! \brief update the state from the pressure
  //! \param [in,out] u   The state to update
  //! \param [in]     eos The state to update
  //============================================================================
  template <typename E>
  static void update_state_from_pressure( state_ref_t & u, const E & eos )
  {
    using math::get;

    // access independant or derived quantities 
    auto m = mass( u );
    auto v = volume( u );
    auto p = pressure( u );
      
    assert( m > 0  );
    assert( v > 0  );
    assert( p > 0  );

    // can use aliases for clarity
    auto & d  = get<variables::index::density>( u );
    auto & ie = get<variables::index::internal_energy>( u );
    auto & t  = get<variables::index::temperature>( u );
    auto & ss = get<variables::index::sound_speed>( u );

    // explicitly set the individual elements
    d  = m / v;
    ie = eos.compute_internal_energy_dp( d, p );
    t  = eos.compute_temperature_de( d, ie );
    ss = eos.compute_sound_speed_de( d, ie );
  }


  //============================================================================
  //! \brief update the state from the energy
  //! \param [in,out] u   The state to update
  //! \param [in]     eos The state to update
  //============================================================================
  template <typename E>
  static void update_state_from_energy( state_ref_t & u, const E & eos )
  {
    using math::get;

    // access independant or derived quantities 
    auto d = density( u );
    auto ie = internal_energy( u );
      
    assert( d > 0  );
    assert( ie > 0  );

    // can use aliases for clarity
    auto & p  = get<variables::index::pressure>( u );
    auto & t  = get<variables::index::temperature>( u );
    auto & ss = get<variables::index::sound_speed>( u );

    // explicitly set the individual elements
    p  = eos.compute_pressure_de( d, ie );
    t  = eos.compute_temperature_de( d, ie );
    ss = eos.compute_sound_speed_de( d, ie );
  }


  //============================================================================
  //! \brief apply an update from conservative fluxes
  //! \param [in,out] u   The state to update
  //! \param [in]     du  The conservative change in state
  //============================================================================
  static void update_state_from_flux( state_ref_t & u, const flux_data_t & du )
  {
    using math::get;
    using math::abs;

    // get the conserved quatities
    auto   m   = mass( u );
    auto   svol= 1 / density( u );
    auto & vel = get<variables::index::velocity>( u );
    auto   et  = total_energy( u );

    // access independant or derived quantities 
    svol += get<equations::index::mass    >( du ) / m;
    vel  += get<equations::index::momentum>( du ) / m;
    et   += get<equations::index::energy  >( du ) / m;

    // get aliases for clarity
    auto & d  = get<variables::index::density>( u );
    auto & ie = get<variables::index::internal_energy>( u );
    auto & vol= get<variables::index::volume>( u );

    // recompute solution quantities
    d = 1 / svol;
    ie = et - 0.5 * abs( vel );
    vol  = m * svol;

    assert( d > 0  );
    assert( ie > 0  );

  }


  //============================================================================
  //! \brief return the volumetric rate of change from the residuals
  //! \param [in]     du  The conservative change in state
  //============================================================================
  static auto volumetric_rate_of_change( const flux_data_t & dudt )
  {
    using math::get;
    return get<equations::index::mass>( dudt );
  }

};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
