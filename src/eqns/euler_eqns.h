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


  //! \brief the variables in the primitive state
  enum equations : size_t { 
    mass = 0, 
    momentum, 
    energy,
    number
  };

  //! \brief the number of equations
  static constexpr size_t num_equations(void)
  {  return equations::number; }
  
  //! \brief return the equation name for a specific component
  //! \param [in] i the component
  //! \return a string with the variable name
  static constexpr std::string eqn_name(auto i) 
  {
    switch ( i ) {
    case equations::mass:
      return "mass";
    case equations::momentum:
      return "momentum";
    case equations::energy:
      return "energy";
    default:
      assert(false && "invalid equation id");
    }
  }
  
  //! \brief return the number of components for the variable
  //! \param [in] i the component
  //! \return the number a components for this variable (1 for scalar, N for vector)    
  static constexpr size_t eqn_components(auto i) 
  {
    switch ( i ) {
    case equations::mass:
      return 1;
    case equations::momentum:
      return dimensions;
    case equations::energy:
      return 1;
    default:
      assert(false && "invalid equation id");
    }
  }
  
  //============================================================================
  //! \brief The primitive state class
  //============================================================================
  class state_t {
  public:

    //! \brief the variables in the primitive state
    enum variables : size_t
    { 
      density = 0, 
      velocity, 
      pressure,
      internal_energy,
      temperature,
      sound_speed,
      number
    };

    //! \brief the number of variables
    static constexpr size_t num_variables(void)
    {  return variables::number; }
    
    //! \brief return the variable name for a specific component
    //! \param [in] i the component
    //! \return a string with the variable name
    static constexpr std::string var_name(auto i) 
    {
      switch ( i ) {
      case variables::density:
        return "density";
      case variables::velocity:
        return "velocity";
      case variables::pressure:
        return "pressure";
      case variables::internal_energy:
        return "internal_energy";
      case variables::temperature:
        return "temperature";
      case variables::sound_speed:
        return "sound_speed";
      default:
        assert(false && "invalid state variable id");
      }
    }
    
    //! \brief return the number of components for the variable
    //! \param [in] i the component
    //! \return the number a components for this variable (1 for scalar, N for vector)    
    static constexpr size_t var_components(auto i) 
    {
      switch ( i ) {
      case variables::density:
        return 1;
      case variables::velocity:
        return eqns_t::dimensions;
      case variables::pressure:
        return 1;
      case variables::internal_energy:
        return 1;
      case variables::temperature:
        return 1;
      case variables::sound_speed:
        return 1;
      default:
        assert(false && "invalid state variable id");
      }
    }

  };



  //! \brief compute the flux in the normal direction
  //! \param [in] normal The normal direction
  //! \return the flux alligned with the normal direction
  static flux_data_t flux( const state_data_t& u, const vector_t& normal )
  {
      
    auto rho = u.template get<state_t::variables::density >();
    auto vel = u.template get<state_t::variables::velocity>();
    auto p   = u.template get<state_t::variables::pressure>();
    auto ie  = u.template get<state_t::variables::internal_energy>();
      
    assert( rho > 0  );

    using math::dot_product;
    auto v_dot_n = dot_product(vel, normal);
    auto pn = p*normal;
    auto et = ie + 0.5 * dot_product(vel, vel);
      
    auto mass_flux = rho * v_dot_n;     
    auto mom_flux  = mass_flux * vel + pn;
    auto ener_flux = (rho*et + p) * v_dot_n;

    return flux_data_t{mass_flux, mom_flux, ener_flux};
  }


};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
