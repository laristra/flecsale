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

//! user includes
#include "ale/common/types.h"
#include "ale/utils/vector.h"

namespace ale {
namespace eqns {

////////////////////////////////////////////////////////////////////////////////
//! \brief Specialization of the euler equations
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
class euler_eqns_t {


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


  //============================================================================
  //! \brief The conserved state class
  //============================================================================
  class conserved_state_t {
  public:

    //! \brief the variables in the conserved state
    enum variables
    { 
      density = 0, 
      momentum, 
      total_energy,
      num_variables };
    
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

    //! \brief Default  constructor
    conserved_state_t() {}

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    conserved_state_t(std::initializer_list<real_t> list) {
      assert(list.size() == variables::num_variables && "dimension size mismatch");
      std::copy(list.begin(), list.end(), data_.begin());
    }

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    template <typename... A> 
    conserved_state_t(A... args) {
      static_assert( (sizeof...(A) == variables::num_variables),
                     "dimension size mismatch" );
      data_ = {args...};
    }


    //! \brief brackets operator
    //! \param [in] i the index
    //! \return the data
    auto &operator[](size_t i)       { return data_[i]; }
    auto  operator[](size_t i) const { return data_[i]; }


    //! \brief Default  constructor (disabled)
    //conserved_state_t() = delete;

    //! \brief Copy constructor (disabled)
    conserved_state_t(const conserved_state_t &) = delete;

    //! \brief Assignment operator (disabled)
    conserved_state_t &operator=(const conserved_state_t &) = delete;

  private:

    //! \brief the data storage
    std::tupel<real_t,vector_t,real_t> data_;

  };
  

  //============================================================================
  //! \brief The primitive state class
  //============================================================================
  class primitive_state_t {
  public:

    //! \brief the variables in the primitive state
    enum variables
    { 
      density = 0, 
      velocity, 
      pressure,
      num_variables 
    };
    
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


    //! \brief Default  constructor
    primitive_state_t() {}

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    primitive_state_t(std::initializer_list<real_t> list) {
      assert(list.size() == variables::num_variables && "dimension size mismatch");
      std::copy(list.begin(), list.end(), data_.begin());
    }

    //! \brief Constructor with initializer list
    //! \param[in] list the initializer list of values
    template <typename... A> 
    primitive_state_t(A... args) {
      static_assert( (sizeof...(A) == variables::num_variables),
                     "dimension size mismatch" );
      data_ = {args...};
    }

    //! \brief brackets operator
    //! \param [in] i the index
    //! \return the data
    template<size_t i>
    auto & get() 
    { return get<i>(data_); }

    template<size_t i>
    auto get() const 
    { return get<i>(data_); }

    //! \brief Default  constructor (disabled)    
    primitive_state_t() = delete;

    //! \brief Copy constructor (disabled)
    primitive_state_t(const primitive_state_t &) = delete;

    //! \brief Assignment operator (disabled)
    primitive_state_t &operator=(const primitive_state_t &) = delete;


  private:

    //! \brief the data storage
    std::tupel<real_t,vector_t,real_t> data_;

  };


  //============================================================================
  //! \brief default constructor
  //============================================================================
  euler_eqns_t() {}


  //============================================================================
  //! \brief conserved to primitive conversion
  //! \param [in]  u the conserved state
  //! \param [out] w the primitive state
  //! \param [in]  get_pressure A function to return the pressure.
  //============================================================================
  void conserved_to_primitive( const conserved_state_t &u, 
                               primitive_state_t &w, 
                               function_t get_pressure ) const 
  {

    const auto & rho = u[conserved_state_t::variables::density];
    const auto & mom = u[conserved_state_t::variables::momentum];
    const auto & et  = u[conserved_state_t::variables::total_energy];

    assert( rho > 0  );

    w[primitive_state_t::variables::density] = rho;
    auto & vel = w[primitive_state_t::variables::velocity];
    auto & p   = w[primitive_state_t::variables::pressure];
 
   
    vel = mom / rho;

    auto e = et/rho - 0.5 * dot(vel,vel);
    p = get_pressure( rho, e );

    
  }

  //============================================================================
  //! \brief primitive to conserved conversion
  //! \param [in]  w the primitive state
  //! \param [out] u the conserved state
  //! \param [in]  get_pressure A function to return the internal energy.
  //============================================================================
  void primitive_to_conserved( const primitive_state_t &w, 
                               conserved_state_t &u,                                
                               function_t get_internal_energy ) const
  {
    
    const auto & rho = w[primitive_state_t::variables::density];
    const auto & vel = w[primitive_state_t::variables::velocity];
    const auto & p   = w[primitive_state_t::variables::pressure];

    assert( rho > 0  );

    u[conserved_state_t::variables::density] = rho;
    auto & mom = u[conserved_state_t::variables::momentum];
    auto & et  = u[conserved_state_t::variables::total_energy];

    mom = rho * vel;

    auto e = get_internal_energy( rho, p );
    et = e + 0.5 * dot(vel, vel);

  }

  //============================================================================
  // Deleted functions
  //============================================================================


  //! \brief disallow copy constructor
  euler_eqns_t(const euler_eqns_t&) = delete;
  
  //! \brief Disallow copying
  euler_eqns_t& operator=(const euler_eqns_t&) = delete;



};

} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
