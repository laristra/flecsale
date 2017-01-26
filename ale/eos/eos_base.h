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
class eos_base_t {


public:


  //============================================================================
  // Typedefs
  //============================================================================

  // \brief the real type
  using real_t = T;


  //============================================================================
  // Public member functions that are part of the common interface
  //============================================================================

  //! \brief return the density
  //! \return the density
  virtual real_t get_ref_density( void ) const = 0;

  //! \brief return the internal energy
  //! \return the internal energy
  virtual real_t get_ref_internal_energy( void ) const = 0;

  //! \brief return the reference temperature
  //! \return the reference temperature
  virtual real_t get_ref_temperature( void ) const = 0;

  //! \brief return the reference pressure
  //! \return the reference pressure
  virtual real_t get_ref_pressure( void ) const = 0;

  //! \brief set the reference state via density and energy
  //! \param[in] density the density to set
  //! \param[in] internal_energy the internal energy to set
  virtual void set_ref_state_de( real_t density, real_t internal_energy ) = 0;

  //! \brief set the reference state via density and temperature
  //! \param[in] density the density to set
  //! \param[in] temperature the temperature to set
  virtual void set_ref_state_dt( real_t density, real_t temperature ) = 0;

  //! \brief set the reference state via density and pressure
  //! \param[in] density the density to set
  //! \param[in] pressure the pressure to set
  virtual void set_ref_state_dp( real_t density, real_t pressure ) = 0;

  //! \brief set the reference state via pressure and temperature
  //! \param[in] pressure the pressure to set
  //! \param[in] temperature the temperature to set
  virtual void set_ref_state_tp( real_t pressure, real_t temperature ) = 0;

  //! \brief compute the internal energy
  //!
  //! \param[in] density the density
  //! \param[in] pressure the pressure
  //! \return the internal energy
  virtual real_t compute_internal_energy_dp( 
    real_t density, 
    real_t pressure 
  ) const = 0;

  //! \brief compute the pressure
  //!
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the pressure
  virtual real_t compute_pressure_de( 
    real_t density, 
    real_t internal_energy 
  ) const = 0;

  //! \brief comput the sound speed
  //!
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the sound speed
  virtual real_t compute_sound_speed_de( 
    real_t density, 
    real_t internal_energy 
  ) const = 0;

  //! \brief comput the temperature
  //!
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the temperature
  virtual real_t compute_temperature_de( 
    real_t density, 
    real_t internal_energy 
  ) const = 0;

  //! \brief Return the gas constant or an effective one. 
  //!
  //! \param[in] density the density
  //! \param[in] pressure the pressure
  //! \return the effective gamma
  virtual real_t compute_gamma_dp( 
    real_t density, 
    real_t pressure 
  ) const = 0;

  //! \brief Return the gas constant or an effective one. 
  //!
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the effective gamma
  virtual real_t compute_gamma_de( 
    real_t density, 
    real_t internal_energy 
  ) const = 0;

};

} // namespace
} // namespace
