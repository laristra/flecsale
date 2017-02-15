/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file lagrange_eqns.h
/// 
/// \brief The desrciption of the euler equations in a lagrangian reference 
///        frame.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "math/tuple.h"
#include "math/general.h"
#include "math/vector.h"
#include "utils/const_string.h"

namespace ale {
namespace eqns {

////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the euler equations in a lagrangian reference 
//!        frame.
//! \tparam T The real type.
//! \tparam N The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T, size_t N>
struct lagrange_eqns_t {


public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief The size type.
  using size_t = std::size_t;

  // \brief The real type.
  using real_t = T;

  //! \brief The vector type.
  using vector_t = math::vector<real_t,N>;

  //! The number of dimensions.
  static constexpr size_t dimensions = N;

  //! A minimum sound speed.
  static constexpr real_t min_sound_speed = 1.e-6;

  //============================================================================
  //! \brief The equations struct.
  //============================================================================
  struct equations {

    //! \brief the variables in the primitive state
    enum index : size_t { 
      mass = 0, 
      momentum, 
      energy = 1 + N,
      total
    };

    //! \brief  The type for holding the state data (mass, momentum, and energy).
    //! 
    //! data_t is a std::tuple with a real_t for mass, a vector_t for momentum
    //! and a real_t for energy.  This needs to correspond to index
    //! or there may be problems
    using data_t = math::array<real_t, index::total>;

    //! \brief the number of equations
    static constexpr size_t number(void)
    {  return index::total; }

  };

  //============================================================================
  //! \brief The variables struct.
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

    //! \brief the number of variables
    static constexpr std::array< utils::const_string, number() > names = 
      { 
        "volume", 
        "mass",
        "velocity", 
        "pressure", 
        "density",
        "internal_energy", 
        "temperature", 
        "sound_speed"
      };
    
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


  //============================================================================
  //! \brief Accessors for various quantities.
  //!
  //! Wrapper functions are only used for accessing either independant or 
  //! derived quantities.  Modifictions of variables should access
  //! the element of the state_data_t explicitily, because it should remain 
  //! clear exactly what variable is being set.  Does 
  //!   set_density(val) 
  //! set density explicitly? Or does it really set temperature and pressure?
  //! But it is clear thatstd::forward<U>(u)
  //!   get<variables::index::density> = val
  //! is assigning a new value to density.
  //!
  //! \param [in] u the state
  //! \return the quantity of interest
  //============================================================================
  template< typename U >
  inline static decltype(auto) volume( U && u )
  { 
    return math::get<variables::index::volume>( std::forward<U>(u) ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) mass( U && u )
  { 
    return math::get<variables::index::mass>( u ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) velocity( U && u )
  { 
    return math::get<variables::index::velocity>( u ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) pressure( U && u )
  { 
    return math::get<variables::index::pressure>( u ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) density( U && u )
  { 
    return math::get<variables::index::density>( u ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) internal_energy( U && u )
  { 
    return math::get<variables::index::internal_energy>( u ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) sound_speed( U && u )
  { 
    return math::get<variables::index::sound_speed>( u ); 
  }

  //! \copydoc volume
  template< typename U >
  inline static decltype(auto) temperature( U && u )
  { 
    return math::get<variables::index::temperature>( u ); 
  }

  //! \copydoc volume
  template<typename U>
  static auto total_energy( U && u )
  { 
    using math::abs;
    auto ie = internal_energy( std::forward<U>(u) );
    auto & vel = velocity( std::forward<U>(u) );
    return ie + dot_product( vel, vel ) / 2;
  }

  //! \copydoc volume
  template<typename U>
  static auto impedance( U && u )
  { 
    return density(std::forward<U>(u)) * sound_speed(std::forward<U>(u));
  }

  //! \copydoc volume
  template <typename U, typename E>
  static auto impedance_multiplier( U && u, const E & eos )
  { 
    // (gamma + 1) / 2
    auto gamma = eos.compute_gamma_de( 
      density(std::forward<U>(u)), internal_energy(std::forward<U>(u)) );
    return 0.5 * ( 1 + gamma );
  }

  //============================================================================
  //! \brief Update the state from the pressure.
  //! \param [in,out] u   The state to update.
  //! \param [in]     eos The equation of state object to apply.
  //! \tparam E  The EOS object type.
  //============================================================================
  template <typename U, typename E>
  static void update_state_from_pressure( U && u, const E & eos )
  {
    using math::get;

    // access independant or derived quantities 
    auto m = mass( std::forward<U>(u) );
    auto v = volume( std::forward<U>(u) );
    auto p = pressure( std::forward<U>(u) );
      
    assert( m > 0  );
    assert( v > 0  );
    assert( p > 0  );

    // can use aliases for clarity
    auto & d  = density( std::forward<U>(u) );
    auto & ie = internal_energy( std::forward<U>(u) );
    auto & t  = temperature( std::forward<U>(u) );
    auto & ss = sound_speed( std::forward<U>(u) );

    // explicitly set the individual elements
    d  = m / v;
    ie = eos.compute_internal_energy_dp( d, p );
    t  = eos.compute_temperature_de( d, ie );
    ss = eos.compute_sound_speed_de( d, ie );
    ss = std::max( ss, min_sound_speed );
  }


  //============================================================================
  //! \brief Update the state from the energy.
  //! \param [in,out] u   The state to update.
  //! \param [in]     eos The equation of state object to apply.
  //! \tparam E  The EOS object type.
  //============================================================================
  template <typename U, typename E>
  static void update_state_from_energy( U && u, const E & eos )
  {
    using math::get;

    // access independant or derived quantities 
    auto d = density( std::forward<U>(u) );
    auto ie = internal_energy( std::forward<U>(u) );
      
    assert( d > 0  );
    assert( ie > 0  );

    // can use aliases for clarity
    auto & p  = pressure( std::forward<U>(u) );
    auto & t  = temperature( std::forward<U>(u) );
    auto & ss = sound_speed( std::forward<U>(u) );

    // explicitly set the individual elements
    p  = eos.compute_pressure_de( d, ie );
    t  = eos.compute_temperature_de( d, ie );
    ss = eos.compute_sound_speed_de( d, ie );
    ss = std::max( ss, min_sound_speed );
  }


  //============================================================================
  //! \brief Apply an update from conservative fluxes.
  //! \param [in,out] u   The state to update.
  //! \param [in]     du  The conservative change in state.
  //============================================================================
  template< typename U, typename F >
  static void update_state_from_flux( 
    U && u, F && du, 
    const real_t & fact = 1.0
  ) {
    using math::get;

    // get the conserved quatities
    auto   m   = mass( std::forward<U>(u) );
    auto & vel = velocity( std::forward<U>(u) );
    auto & ie  = internal_energy( std::forward<U>(u) );

    // get the changes
    auto inv_mass = 1 / m;    
    
    // update total energy
    auto et = total_energy(std::forward<U>(u)) + 
              fact * inv_mass * std::forward<F>(du)[equations::index::energy];

    // update velocity
    for ( int i=0; i<N; ++i )
      vel[i] += fact * inv_mass * 
                std::forward<F>(du)[equations::index::momentum + i];

    // compute new internal
    ie = et - dot_product( vel, vel ) / 2;

    if ( internal_energy(u) < 0 )
      raise_runtime_error( 
        "Negative internal energy encountered, " << internal_energy(u) << "." 
        << std::endl << "Current state = " << u << "."
      );

  }


  //============================================================================
  //! \brief Update the volume, and consequentially the density assuming mass
  //!        remains constant..
  //! \param [in,out] u   The state to update.
  //! \param [in]     du  The conservative change in state.
  //============================================================================
  template< typename U >
  static void update_volume( U && u, real_t new_vol )
  {
    using math::get;

    // recompute solution quantities
    volume(std::forward<U>(u)) = new_vol;
    density(std::forward<U>(u)) = mass(std::forward<U>(u)) / new_vol;

    if ( density(std::forward<U>(u)) < 0 )
      raise_runtime_error( 
        "Negative density encountered, " << density(u) << "." 
        << std::endl << "Current state = " << u << "."
      );

  }


  //============================================================================
  //! \brief Update the volume, and consequentially the density assuming mass
  //!        remains constant..
  //============================================================================
  template< typename U >
  static void compute_update( 
    const vector_t & u,      //!< [in] the point velocity
    const vector_t & force,  //!< [in] the force at the point
    const vector_t & n,      //!< [in] the normal direction to apply the force
    U && dudt                //!< [in,out] the residual to update
  ) {

      // add contribution
      for ( int d=0; d<N; ++d ) {
        std::forward<U>(dudt)[equations::index::mass] += n[d] * u[d];
        std::forward<U>(dudt)[equations::index::energy] -= force[d] * u[d];
        std::forward<U>(dudt)[equations::index::momentum + d] -= force[d];
      }

  }

  //============================================================================
  //! \brief Return the volumetric rate of change from the residuals.
  //! \param [in]     du  The conservative change in state.
  //============================================================================
  template< typename U >
  static auto volumetric_rate_of_change( U && dudt )
  {
    return std::forward<U>(dudt)[equations::index::mass];
  }

};



////////////////////////////////////////////////////////////////////////////////
// a minimum sound speed
// note: this has to be here
////////////////////////////////////////////////////////////////////////////////
template<typename T, size_t N>
constexpr T lagrange_eqns_t<T,N>::min_sound_speed;

} // namespace
} // namespace

