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

  static auto pressure( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::pressure>( u ); 
  }

  static auto internal_energy( const state_data_t & u )
  { 
    using math::get;
    return get<variables::index::pressure>( u ); 
  }

  static auto total_energy( const state_data_t & u )
  { 
    using math::get;
    using math::dot_product;
    auto ie = internal_energy( u );
    auto vel = velocity( u );
    return ie + 0.5 * dot_product(vel, vel);
  }

  //============================================================================
  //! \brief compute the flux in the normal direction
  //! \param [in] normal The normal direction
  //! \return the flux alligned with the normal direction
  //============================================================================
  static auto flux( const auto & u, const auto & normal )
  {

    using math::get;
    using math::dot_product;

    // these may be independant or derived quantities
    auto rho = density ( u );
    auto vel = velocity( u );
    auto p   = pressure( u );
    auto et  = total_energy( u );
      
    assert( rho > 0  );

    auto v_dot_n = dot_product(vel, normal);
    auto pn = p*normal;
      
    // explicitly set the individual elements, and it is clear what is
    // being set by using static indexing instead of wrapper functions
    flux_data_t f;
    auto mass_flux = rho * v_dot_n;
    get<equations::index::mass    >( f ) = mass_flux;     
    get<equations::index::momentum>( f ) = mass_flux * vel + pn;
    get<equations::index::energy  >( f ) = (rho*et + p) * v_dot_n;

    return f;
  }


  //============================================================================
  //! \brief update the state from the pressure
  //! \param [in,out] u   The state to update
  //! \param [in]     eos The state to update
  //============================================================================
  static void update_state_from_pressure( auto & u, 
                                          const auto& eos )
  {
    using math::get;

    // access independant or derived quantities 
    auto d = density( u );
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


};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
