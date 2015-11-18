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
#include "ale/common/types.h"
#include "ale/utils/vector.h"

namespace ale {
namespace eqns {


////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the euler equations
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
struct euler_eqns_t {


public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the real precision type
  using real_t = common::real_t;

  //! \brief the size type
  using size_t = common::size_t;

  //! \brief the vector_type
  using vector_t = utils::vector_t<real_t,N>;

  //! \brief a function template for the eos closure
  using function_t = std::function<real_t(real_t, real_t)>;

  //! the number of dimensions
  static constexpr size_t dimensions = N;
  //! the number of equations to solve for
  static constexpr size_t equations = 2 + N;

  //! this equation set has a primitive state
  static constexpr bool has_primitive_state = true;

  //! \brief an alias to the eos type
  using eqns_t = euler_eqns_t;


  //! \brief  the type for holding the state data
  using state_data_t = std::tuple<real_t,vector_t,real_t>;


private:

  //============================================================================
  //! \brief The base state class
  //============================================================================
  class base_state_t {
  public:


    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    template <typename... Args>
    base_state_t(Args&&... args) : data_(std::forward<Args>(args)...) 
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


    //! \brief Output operator
    //! \param[in,out] os  The ostream to dump output to.
    //! \param[in]     rhs The vector_t on the right hand side of the operator.
    //! \return A reference to the current ostream.
    friend auto & operator<<(std::ostream& os, const base_state_t& a)
    {
      os << "{";
      os << " "  << a.template get<0>();
      os << ", " << a.template get<1>();
      os << ", " << a.template get<2>();
      os << " }";
      return os;
    }



  protected:

    //! \brief the data storage
    state_data_t data_;

  };

public:

  //! forward decares
  class conserved_state_t;
  class primitive_state_t;
  

  //============================================================================
  //! \brief The conserved state class
  //============================================================================
  class conserved_state_t : public base_state_t {
  public:

    //! \brief the variables in the conserved state
    enum variables : size_t
    { 
      density = 0, 
      momentum, 
      total_energy,
      number };

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
      case variables::momentum:
        return "momentum";
      case variables::total_energy:
        return "total_energy";
      default:
        assert(false && "invalid conserved variable id");
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
      case variables::momentum:
        return eqns_t::dimensions;
      case variables::total_energy:
        return 1;
      default:
        assert(false && "invalid conserved variable id");
      }
    }

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    template <typename... Args,
              typename = typename std::enable_if<sizeof...(Args) == num_variables()>::type>
    conserved_state_t(Args&&... args) : base_state_t(std::forward<Args>(args)...) 
    { }

    //! \brief Constructor with state_data_
    //! \param[in] list the initializer list of values
    conserved_state_t(const state_data_t & other) : base_state_t( other ) 
    {  }
    
    //! \brief conserved to primitive conversion
    //! \param [in]  get_pressure A function to return the pressure.
    //! \return the primitive state
    primitive_state_t to_primitive( function_t get_pressure ) const 
    {
      auto rho = base_state_t::template get<variables::density>();
      auto mom = base_state_t::template get<variables::momentum>();
      auto et  = base_state_t::template get<variables::total_energy>();
      
      assert( rho > 0  );
      
      auto vel = mom / rho;
      auto e   = et/rho - 0.5 * utils::dot(vel,vel);
      auto p   = get_pressure( rho, e );
      
      return primitive_state_t{rho,vel,p};
    }

    //! \brief Equivalence operator
    //! \param[in] lhs The quantity on the rhs.
    //! \param[in] rhs The quantity on the rhs.
    //! \return true if equality.
    friend bool operator==(const conserved_state_t& lhs, const conserved_state_t& rhs)
    { 
      return (lhs.data_ == rhs.data_); 
    }


  };
  

  //============================================================================
  //! \brief The primitive state class
  //============================================================================
  class primitive_state_t : public base_state_t {
  public:

    //! \brief the variables in the primitive state
    enum variables : size_t
    { 
      density = 0, 
      velocity, 
      pressure,
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
      default:
        assert(false && "invalid conserved variable id");
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
      default:
        assert(false && "invalid conserved variable id");
      }
    }

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    template <typename... Args,
              typename = typename std::enable_if<sizeof...(Args) == num_variables()>::type>
    primitive_state_t(Args&&... args) : base_state_t(std::forward<Args>(args)...) 
    {  }

    //! \brief Constructor with state_data_
    //! \param[in] list the initializer list of values
    primitive_state_t(const state_data_t & other) : base_state_t( other ) 
    {  }


    //! \brief primitive to conserved conversion
    //! \param [in]  get_pressure A function to return the internal energy.
    //! \return the conserved state
    conserved_state_t to_conserved( function_t get_internal_energy ) const
    {
      
      auto rho = base_state_t::template get<variables::density >();
      auto vel = base_state_t::template get<variables::velocity>();
      auto p   = base_state_t::template get<variables::pressure>();
      
      assert( rho > 0  );
      
      auto mom = rho * vel;     
      auto e  = get_internal_energy( rho, p );
      auto et = e + 0.5 * utils::dot(vel, vel);
      
      return conserved_state_t{rho, mom, et};
    }
    

    //! \brief compute the flux in the normal direction
    //! \param [in] normal The normal direction
    //! \return the flux alligned with the normal direction
    state_data_t flux( const vector_t normal ) const
    {
      
      auto rho = base_state_t::template get<variables::density >();
      auto vel = base_state_t::template get<variables::velocity>();
      auto p   = base_state_t::template get<variables::pressure>();
      
      assert( rho > 0  );

      auto v_dot_n = utils::dot(vel, normal);
      
      auto mass_flux = rho * v_dot_n;     
      auto mom_flux  = mass_flux * vel;
      auto ener_flux = mom_flux + p;
      
      return std::make_tuple(mass_flux, mom_flux, ener_flux);
    }

    //! \brief Equivalence operator
    //! \param[in] lhs The quantity on the rhs.
    //! \param[in] rhs The quantity on the rhs.
    //! \return true if equality.
    friend bool operator==(const primitive_state_t& lhs, const primitive_state_t& rhs)
    { 
      return (lhs.data_ == rhs.data_); 
    }

  };

};


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
