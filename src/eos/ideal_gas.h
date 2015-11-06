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
 * \file ideal_gas.h
 * 
 * \brief Ideal gas equation of state implementation.
 *
 ******************************************************************************/
#pragma once

//! system includes

//! user includes
#include "ale/common/types.h"
#include "ale/utils/zip.h"

namespace ale {
namespace eos {

////////////////////////////////////////////////////////////////////////////////
//! \brief Ideal gas specialization of the equation of state
////////////////////////////////////////////////////////////////////////////////
class ideal_gas_t {


public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the real precision type
  using real_t = common::real_t;

  //============================================================================
  // Constructors / Destructors
  //============================================================================


  //! \brief default constructor
  ideal_gas_t() : gamma_(2.0), density_(1.0), specific_heat_v_(1.0), 
                  internal_energy_(1.0)
  {}

  //! \brief the destructor
  ~ideal_gas_t() {};

  //! \brief disallow copy constructor
  ideal_gas_t(const ideal_gas_t&) = delete;
  
  //! \brief Disallow copying
  ideal_gas_t& operator=(const ideal_gas_t&) = delete;

  //============================================================================
  // Public member functions that are special for this class 
  //============================================================================

  //! \brief return the specific heat at constant volume
  //! \return the specific heat
  real_t get_specific_heat_v( void ) const 
  { return specific_heat_v_; }

  //! \brief set the specific heat of constant volume
  //! \param[in] specific_heat_v the specific heat of constant volume to set
  void set_specific_heat_v( real_t specific_heat_v )
  { specific_heat_v_ = specific_heat_v; }


  //! \brief return the specific heat ratio
  //! \return the specific heat ratio
  real_t get_gamma( void ) const 
  { return gamma_; }

  //! \brief set the specific heat ratio
  //! \param[in] specific_heat_v the specific heat ratio to set
  void set_gamma( real_t gamma )
  { gamma_ = gamma; }


  //============================================================================
  // Public member functions that are part of the common interface
  //============================================================================

  //! \brief Main driver to setup the eos using input data
  //! \param[in] input the input object with the data
  //virtual void initialize( const common::NodeType &input ) = 0;

  //! \brief Main driver to setup the eos object
  void setup( void ) {};


  //! \brief return the density
  //! \return the density
  real_t get_ref_density( void ) const 
  { return density_; }

  //! \brief set the density
  //! \param[in] density the density to set
  void set_ref_density( real_t density )
  { density_ = density; }

  //! \brief return the internal energy
  //! \return the internal energy
  real_t get_ref_internal_energy( void ) const 
  { return internal_energy_; }

  //! \brief set the internal energy
  //! \param[in] internal_energy the internal energy to set
  void set_ref_internal_energy( real_t internal_energy )
  { internal_energy_ = internal_energy; }


  //! \brief compute the internal energy
  //!
  //! This version is less efficient since there is one virtual call for
  //! every internal energy evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] pressure the pressure
  //! \return the internal energy
  real_t compute_internal_energy( real_t density, 
                                  real_t pressure ) const
  { return pressure / ( density*(gamma_-1.0) ); }


  //! \brief compute the internal energy
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of internal energies.
  //! 
  //! \param[in]  density the density
  //! \param[in]  pressure the pressure
  //! \param[out] internal_energy the internal energy
  void compute_internal_energy( const auto & density, 
                                const auto & pressure,
                                auto & internal_energy ) const
  {
    auto tup = utils::zip( density, pressure, internal_energy );
    for ( auto it : tup ) {
      auto & d = it.template get<0>();
      auto & p = it.template get<1>();
      auto & e = it.template get<2>();
      e = compute_internal_energy( d, p );
    }
  }


  //! \brief compute the pressure
  //!
  //! This version is less efficient since there is one virtual call for
  //! every internal pressure evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the pressure
  real_t compute_pressure( real_t density, 
                           real_t internal_energy ) const
  { return (gamma_-1.0) * density * internal_energy; }

  //! \brief compute the pressure
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of pressures.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \param[out] pressure the pressure
  void compute_pressure( const auto & density, 
                         const auto & internal_energy,
                         auto & pressure ) const
  {
    auto tup = utils::zip( density, internal_energy, pressure );
    for ( auto it : tup ) {
      auto & d = it.template get<0>();
      auto & e = it.template get<1>();
      auto & p = it.template get<2>();
      p = compute_pressure( d, e );
    }
  }



  //! \brief comput the sound speed
  //!
  //! This version is less efficient since there is one virtual call for
  //! every sound speed evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the sound speed
  real_t compute_sound_speed( real_t density, 
                              real_t internal_energy ) const
  {
    assert( internal_energy>0.0 );
    return std::sqrt( gamma_ * (gamma_-1.0) * internal_energy );
  }

  //! \brief comput the sound speed
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of sound speeds.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \param[out] sound_speed the sound speed
  void compute_sound_speed( const auto & density, 
                            const auto & internal_energy,
                            auto & sound_speed ) const
  {
    auto tup = utils::zip( density, internal_energy, sound_speed );
    for ( auto it : tup ) {
      auto & d = it.template get<0>();
      auto & e = it.template get<1>();
      auto & s = it.template get<2>();
      s = compute_sound_speed( d, e );
    }
  }


  //! \brief comput the temperature
  //!
  //! This version is less efficient since there is one virtual call for
  //! every temperature evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the temperature
  real_t compute_temperature( real_t density, 
                              real_t internal_energy ) const
  { return internal_energy / specific_heat_v_; }


  //! \brief comput the temperature
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of temperatures.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \param[out] temperature the temperature
  void compute_temperature( const auto & density, 
                            const auto & internal_energy,
                            auto & temperature ) const
  {
    auto tup = utils::zip( density, internal_energy, temperature );
    for ( auto it : tup ) {
      auto & d = it.template get<0>();
      auto & e = it.template get<1>();
      auto & t = it.template get<2>();
      t = compute_temperature( d, e );
    }
  }


  protected:
    
  //===============================================================
  // Member variables
  //===============================================================
  
  //! the ideal gas specific heat ration
  real_t gamma_; 

  //! the specific heat
  real_t specific_heat_v_; 

  //! the reference density
  real_t density_; 

  //! the reference internal energy
  real_t internal_energy_; 


};

} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
