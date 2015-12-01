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

//! system includes
#include <functional>
#include <iostream>
#include <tuple>

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

  //! \brief a function template for the eos closure
  using function_t = std::function<real_t(real_t, real_t)>;

  //! the number of dimensions
  static constexpr size_t dimensions = N;


  //! \brief an alias to the eos type
  using eqns_t = euler_eqns_t;

  //! \brief the vector_type
  using vector_t = math::vector_t<real_t,N>;

  //! \brief  the type for holding the state data
  using state_data_t = math::tuple_t<real_t,vector_t,real_t,real_t,real_t,real_t>;

  //! \brief  the type for holding the state data
  using flux_data_t = math::tuple_t<real_t,vector_t,real_t>;

  //============================================================================
  // \brief The equations struct
  //============================================================================
  struct equations {

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
  
    //! \brief return the equation name for a specific component
    //! \param [in] i the component
    //! \return a string with the variable name
    static constexpr const char * name(auto i) 
    {
      switch ( i ) {
      case index::mass:
        return "mass";
      case index::momentum:
        return "momentum";
      case index::energy:
        return "energy";
      }
    }
  
    //! \brief return the number of components for the variable
    //! \param [in] i the component
    //! \return the number a components for this variable (1 for scalar, N for vector)    
    static constexpr size_t components(auto i) 
    {
      switch ( i ) {
      case index::mass:
        return 1;
      case index::momentum:
        return dimensions;
      case index::energy:
        return 1;
      default:
        assert(false && "invalid equation id");
      }
    }
    
  };

  //============================================================================
  // \brief The variables struct
  //============================================================================
  struct variables {

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
    
    //! \brief return the variable name for a specific component
    //! \param [in] i the component
    //! \return a string with the variable name
    static constexpr const char * name(auto i) 
    {
      switch ( i ) {
      case index::density:
        return "density";
      case index::velocity:
        return "velocity";
      case index::pressure:
        return "pressure";
      case index::internal_energy:
        return "internal_energy";
      case index::temperature:
        return "temperature";
      case index::sound_speed:
        return "sound_speed";
      }
    }
    
    //! \brief return the number of components for the variable
    //! \param [in] i the component
    //! \return the number a components for this variable (1 for scalar, N for vector)    
    static constexpr size_t components(auto i) 
    {
      switch ( i ) {
      case index::density:
        return 1;
      case index::velocity:
        return eqns_t::dimensions;
      case index::pressure:
        return 1;
      case index::internal_energy:
        return 1;
      case index::temperature:
        return 1;
      case index::sound_speed:
        return 1;
      default:
        assert(false && "invalid state variable id");
      }
    }

  };

  //! \brief compute the flux in the normal direction
  //! \param [in] normal The normal direction
  //! \return the flux alligned with the normal direction
  static flux_data_t flux( const auto & u, const auto& normal )
  {
    using math::operator*;
    using math::operator+;

    using std::get;
    auto rho = get<variables::index::density >( u );
    auto vel = get<variables::index::velocity>( u );
    auto p   = get<variables::index::pressure>( u );
    auto ie  = get<variables::index::internal_energy>( u );
      
    assert( rho > 0  );

    using math::dot_product;
    auto v_dot_n = dot_product(vel, normal);
    auto pn = p*normal;
    auto et = ie + 0.5 * dot_product(vel, vel);
      
    auto mass_flux = rho * v_dot_n;     
    auto mom_flux  = mass_flux * vel + pn;
    auto ener_flux = (rho*et + p) * v_dot_n;

    return flux_data_t(mass_flux, mom_flux, ener_flux);
  }


  //! \brief update the state from the pressure
  //! \param [in,out] u   The state to update
  //! \param [in]     eos The state to update
  static void update_state_from_pressure( auto & u, 
                                          const auto& eos )
  {
    using std::get;
    auto &d   = get<variables::index::density >( u );
    auto &p   = get<variables::index::pressure>( u );
    auto &ie  = get<variables::index::internal_energy>( u );
    auto &t   = get<variables::index::temperature>( u );
    auto &ss  = get<variables::index::sound_speed>( u );
      
    assert( d > 0  );
    assert( p > 0  );

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
