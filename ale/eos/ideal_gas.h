/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Ideal gas equation of state implementation.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <cmath>

namespace ale {
namespace eos {

////////////////////////////////////////////////////////////////////////////////
//! \brief Ideal gas specialization of the equation of state
////////////////////////////////////////////////////////////////////////////////
template <typename T>
class ideal_gas_t {


public:


  //============================================================================
  // Typedefs
  //============================================================================

  // \brief the real type
  using real_t = T;

  //============================================================================
  // Constructors / Destructors
  //============================================================================


  //! \brief default constructor
  ideal_gas_t() : gamma_(2.0), specific_heat_v_(1.0), 
                  density_(1.0), internal_energy_(1.0)
  {}

  //! \brief constructor with values constructor
  ideal_gas_t( real_t gamma, real_t cv ) : gamma_(gamma), specific_heat_v_(cv),
                                           density_(1.0), internal_energy_(1.0)
  {
    assert( gamma > 1 );
    assert( cv > 0 );
  }

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
  { 
    assert( specific_heat_v > 0 );
    specific_heat_v_ = specific_heat_v; 
  }


  //! \brief return the specific heat ratio
  //! \return the specific heat ratio
  real_t get_gamma( void ) const 
  { return gamma_; }

  //! \brief set the specific heat ratio
  //! \param[in] specific_heat_v the specific heat ratio to set
  void set_gamma( real_t gamma )
  { 
    assert( gamma > 1 );
    gamma_ = gamma; 
  }


  //============================================================================
  // Public member functions that are part of the common interface
  //============================================================================

  //! \brief Main driver to setup the eos using input data
  //! \param[in] input the input object with the data
  //virtual void initialize( const common::NodeType &input ) = 0;

  //! \brief Main driver to setup the eos object
  //void setup( void ) {};


  //! \brief return the density
  //! \return the density
  real_t get_ref_density( void ) const 
  { return density_; }

  //! \brief return the internal energy
  //! \return the internal energy
  real_t get_ref_internal_energy( void ) const 
  { return internal_energy_; }

  //! \brief return the reference temperature
  //! \return the reference temperature
  real_t get_ref_temperature( void ) const 
  { return compute_temperature_de( density_, internal_energy_ ); }

  //! \brief return the reference pressure
  //! \return the reference pressure
  real_t get_ref_pressure( void ) const 
  { return compute_pressure_de( density_, internal_energy_ ); }


  //! \brief set the reference state via density and energy
  //! \param[in] density the density to set
  //! \param[in] internal_energy the internal energy to set
  void set_ref_state_de( real_t density, real_t internal_energy )
  { 
    assert( density > 0 );
    assert( internal_energy > 0 );

    density_ = density; 
    internal_energy_ = internal_energy;
  }

  //! \brief set the reference state via density and temperature
  //! \param[in] density the density to set
  //! \param[in] temperature the temperature to set
  void set_ref_state_dt( real_t density, real_t temperature )
  { 
    assert( density > 0 );
    assert( temperature > 0 );

    density_ = density; 
    internal_energy_ = specific_heat_v_ * temperature; 
  }

  //! \brief set the reference state via density and pressure
  //! \param[in] density the density to set
  //! \param[in] pressure the pressure to set
  void set_ref_state_dp( real_t density, real_t pressure )
  { 
    assert( density > 0 );
    assert( pressure > 0 );

    density_ = density; 
    internal_energy_ = pressure / ( density_*(gamma_-1.0) );
  }


  //! \brief set the reference state via pressure and temperature
  //! \param[in] pressure the pressure to set
  //! \param[in] temperature the temperature to set
  void set_ref_state_tp( real_t pressure, real_t temperature )
  { 
    assert( pressure > 0 );
    assert( temperature > 0 );

    auto R = (gamma_-1.0) * specific_heat_v_;

    density_ = pressure / ( R * temperature); 
    internal_energy_ = pressure / ( density_*(gamma_-1.0) );
  }


  //! \brief compute the internal energy
  //!
  //! This version is less efficient since there is one virtual call for
  //! every internal energy evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] pressure the pressure
  //! \return the internal energy
  real_t compute_internal_energy_dp( real_t density, 
                                     real_t pressure ) const
  { return pressure / ( density*(gamma_-1.0) ); }



  //! \brief compute the pressure
  //!
  //! This version is less efficient since there is one virtual call for
  //! every internal pressure evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the pressure
  real_t compute_pressure_de( real_t density, 
                              real_t internal_energy ) const
  { return (gamma_-1.0) * density * internal_energy; }

  //! \brief comput the sound speed
  //!
  //! This version is less efficient since there is one virtual call for
  //! every sound speed evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the sound speed
  real_t compute_sound_speed_de( real_t density, 
                                 real_t internal_energy ) const
  {
    assert( internal_energy>0.0 );
    return std::sqrt( gamma_ * (gamma_-1.0) * internal_energy );
  }

  //! \brief comput the temperature
  //!
  //! This version is less efficient since there is one virtual call for
  //! every temperature evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the temperature
  real_t compute_temperature_de( real_t density, 
                                 real_t internal_energy ) const
  { return internal_energy / specific_heat_v_; }


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
