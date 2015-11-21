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
#include "ale/utils/vector.h"

namespace ale {
namespace eqns {


////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the euler equations
////////////////////////////////////////////////////////////////////////////////
template<typename real_t, size_t N>
struct euler_eqns_t {


public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the size type
  using size_t = std::size_t;

  //! \brief a function template for the eos closure
  using function_t = std::function<real_t(real_t, real_t)>;

  //! the number of dimensions
  static constexpr size_t dimensions = N;


  //! \brief an alias to the eos type
  using eqns_t = euler_eqns_t;

  //! \brief the vector_type
  using vector_t = utils::vector_t<real_t,N>;

  //! \brief  the type for holding the state data
  using state_data_t = std::tuple<real_t,vector_t,real_t,real_t,real_t,real_t>;

  //! \brief  the type for holding the state data
  using flux_data_t = std::tuple<real_t,vector_t,real_t>;


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

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    template <typename... Args,
              typename = typename std::enable_if<sizeof...(Args) == num_variables()>::type>
    state_t(Args&&... args) : data_(std::forward<Args>(args)...) 
    { }

    //! \brief brackets operator
    //! \param [in] i the index
    //! \return the data
    template<size_t i>
    auto & get() 
    { return std::get<i>(data_); }

    template<size_t i>
    auto get() const 
    { return std::get<i>(data_); }


    //! \brief compute the flux in the normal direction
    //! \param [in] normal The normal direction
    //! \return the flux alligned with the normal direction
    flux_data_t flux( const vector_t normal ) const
    {
      
      auto rho = this->template get<variables::density >();
      auto vel = this->template get<variables::velocity>();
      auto p   = this->template get<variables::pressure>();
      auto ie  = this->template get<variables::internal_energy>();
      
      assert( rho > 0  );

      auto v_dot_n = utils::dot(vel, normal);
      auto pn = p*normal;
      auto et = ie + 0.5 * utils::dot(vel, vel);
      
      auto mass_flux = rho * v_dot_n;     
      auto mom_flux  = mass_flux * vel + pn;
      auto ener_flux = (rho*et + p) * v_dot_n;

      return flux_data_t{mass_flux, mom_flux, ener_flux};
    }


    //! \brief Output operator
    //! \param[in,out] os  The ostream to dump output to.
    //! \param[in]     rhs The vector_t on the right hand side of the operator.
    //! \return A reference to the current ostream.
    friend auto & operator<<(std::ostream& os, const state_t& a)
    {
      os << "{";
      os << " "  << a.template get<0>();
      os << ", " << a.template get<1>();
      os << ", " << a.template get<2>();
      os << ", " << a.template get<3>();
      os << ", " << a.template get<4>();
      os << ", " << a.template get<5>();
      os << " }";
      return os;
    }


    //! \brief Equivalence operator
    //! \param[in] lhs The quantity on the rhs.
    //! \param[in] rhs The quantity on the rhs.
    //! \return true if equality.
    friend bool operator==(const state_t& lhs, const state_t& rhs)
    { 
      return (lhs.data_ == rhs.data_); 
    }


  protected:

    //! \brief the data storage
    state_data_t data_;

  };


};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
