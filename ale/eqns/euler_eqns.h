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
 * \file euler_eqns.h
 * 
 * \brief The desrciption of the euler equations.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/math/tuple.h"
#include "ale/math/math.h"
#include "ale/math/vector.h"

namespace ale {
namespace eqns {

////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the euler equations
////////////////////////////////////////////////////////////////////////////////
template<typename T, size_t N>
struct euler_eqns_t {


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
      math::tuple<real_t,vector_t,real_t,real_t,real_t,real_t>;

    //! \brief the variables in the primitive state
    enum index : size_t
    { 
      density = 0, 
      velocity, 
      pressure,
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
  static auto density( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::density>( u ); 
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
    using math::get;
    using math::dot_product;
    auto ie = internal_energy( u );
    auto vel = velocity( u );
    return ie + dot_product( vel, vel ) / 2;
  }

  //============================================================================
  //! \brief Compute the fastest moving eigenvalue
  //! \param [in]  u   The solution state
  //! \return the fastest moving wave speed
  //============================================================================
  template <typename V>
  static auto fastest_wavespeed( const state_data_t & u, const V & norm )
  {
    auto a = sound_speed( u );
    auto v = velocity( u );
    auto vn = math::dot_product( v, norm );
    return a + std::abs(vn);

  }

  //============================================================================
  //! \brief Compute the fastest moving eigenvalue
  //! \param [in]  u   The solution state
  //! \return the fastest moving wave speed
  //============================================================================
  template <typename V>
  static auto eigenvalues( const state_data_t & u, const V & norm )
  {
    using math::get;

    auto a = sound_speed( u );
    auto v = velocity( u );
    auto vn = math::dot_product( v, norm );
    flux_data_t eig;
    get<equations::index::mass    >( eig ) = vn-a;
    get<equations::index::momentum>( eig ) = vn;
    get<equations::index::energy  >( eig ) = vn+a;
    return eig;
  }

  //============================================================================
  //! \brief Compute the fastest moving eigenvalue
  //! \param [in]  u   The solution state
  //! \return the fastest moving wave speed
  //============================================================================
  template <typename V>
  static auto minmax_eigenvalues( const state_data_t & u, const V & norm )
  {
    using math::get;

    auto a = sound_speed( u );
    auto v = velocity( u );
    auto vn = math::dot_product( v, norm );
    return std::make_pair( vn-a, vn+a );
  }

  //============================================================================
  //! \brief Computes the change in conserved quantities between two states
  //! \param [in]  ul   The left state
  //! \param [in]  ur   The right state
  //! \return ur - ul
  //============================================================================
  static auto solution_delta( 
    const state_data_t & ul, const state_data_t & ur )
  {
    using math::get;

    // get the conserved quatities
    auto mass_l = density( ul );
    auto mom_l  = mass_l*velocity( ul );
    auto ener_l = mass_l*total_energy( ul );

    auto mass_r = density( ur );
    auto mom_r  = mass_r*velocity( ur );
    auto ener_r = mass_r*total_energy( ur );

    // compute the change
    flux_data_t du;
    get<equations::index::mass    >( du ) = mass_r - mass_l;
    get<equations::index::momentum>( du ) = mom_r - mom_l;
    get<equations::index::energy  >( du ) = ener_r - ener_l;

    return du;

  }

  //============================================================================
  //! \brief compute the flux in the normal direction
  //! \param [in] normal The normal direction
  //! \return the flux alligned with the normal direction
  //============================================================================
  template <typename V>
  static auto flux( const state_data_t & u, const V & norm )
  {

    using math::get;
    using math::dot_product;

    // these may be independant or derived quantities
    auto rho = density ( u );
    auto vel = velocity( u );
    auto p   = pressure( u );
    auto et  = total_energy( u );
      
    assert( rho > 0  );

    auto v_dot_n = dot_product( vel, norm );
    auto pn = p*norm;
      
    // explicitly set the individual elements, and it is clear what is
    // being set by using static indexing instead of wrapper functions
    flux_data_t f;
    auto mass_flux = rho * v_dot_n;
    get<equations::index::mass    >( f ) = mass_flux;     
    get<equations::index::momentum>( f ) = mass_flux * vel + pn;
    get<equations::index::energy  >( f ) = mass_flux * (et + p/rho);

    return f;
  }

  //============================================================================
  //! \brief The flux at the wall
  //!
  //! \param [in] u      the solution state
  //! \param [in] norm  the normal direction
  //! \return the flux
  //============================================================================
  template <typename V>
  static auto wall_flux( const state_data_t & u, const V & norm )
  {
    auto p  = pressure( u );
    auto pn = p*norm;
    return flux_data_t{ 0.0, pn, 0.0 };
  };


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
    auto d = density ( u );
    auto p = pressure( u );
      
    assert( d > 0  );
    assert( p > 0  );


    // can use aliases for clarity
    auto & ie = get<variables::index::internal_energy>( u );
    auto & t  = get<variables::index::temperature>( u );
    auto & ss = get<variables::index::sound_speed>( u );

    // explicitly set the individual elements
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
    auto d  = density ( u );
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
    auto mass = density( u );
    auto mom  = mass*velocity( u );
    auto ener = mass*total_energy( u );

    // access independant or derived quantities 
    mass += get<equations::index::mass    >( du );     
    mom  += get<equations::index::momentum>( du );
    ener += get<equations::index::energy  >( du );

    // get aliases for clarity
    auto & d  = get<variables::index::density>( u );
    auto & v  = get<variables::index::velocity>( u );
    auto & ie = get<variables::index::internal_energy>( u );

    // recompute solution quantities
    d = mass;
    v = mom / mass;
    auto et = ener / mass;
    ie = et - 0.5 * dot_product( v, v );

    assert( d > 0  );
    assert( ie > 0  );

  }

};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
