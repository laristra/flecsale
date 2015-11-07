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
 * \file eos_utils.h
 * 
 * \brief Utilities for creating eos tasks.
 *
 ******************************************************************************/
#pragma once

//! system includes
#include <vector>

//! user includes
#include "ale/utils/zip.h"

namespace ale {
namespace eos {

////////////////////////////////////////////////////////////////////////////////
//! \brief compute the internal energy
//!
//! This version is more efficient since there is only one virtual call to
//! evaluate a list of internal energies.
//! 
//! \param[in]  density the density
//! \param[in]  pressure the pressure
//! \param[out] internal_energy the internal energy
////////////////////////////////////////////////////////////////////////////////
void compute_internal_energy( const auto & eos,
                              const auto & density, 
                              const auto & pressure,
                              auto & internal_energy )
{
  auto tup = utils::zip( density, pressure, internal_energy );
  for ( auto it : tup ) {
    auto & d = it.template get<0>();
    auto & p = it.template get<1>();
    auto & e = it.template get<2>();
    e = eos.compute_internal_energy( d, p );
  }
}

////////////////////////////////////////////////////////////////////////////////
//! \brief compute the pressure
//!
//! This version is more efficient since there is only one virtual call to
//! evaluate a list of pressures.
//! 
//! \param[in] density the density
//! \param[in] internal_energy the internal energy
//! \param[out] pressure the pressure
////////////////////////////////////////////////////////////////////////////////
void compute_pressure( const auto & eos,
                       const auto & density, 
                       const auto & internal_energy,
                       auto & pressure )
{
  auto tup = utils::zip( density, internal_energy, pressure );
  for ( auto it : tup ) {
    auto & d = it.template get<0>();
    auto & e = it.template get<1>();
    auto & p = it.template get<2>();
    p = eos.compute_pressure( d, e );
  }
}




////////////////////////////////////////////////////////////////////////////////
//! \brief comput the sound speed
//!
//! This version is more efficient since there is only one virtual call to
//! evaluate a list of sound speeds.
//! 
//! \param[in] density the density
//! \param[in] internal_energy the internal energy
//! \param[out] sound_speed the sound speed
////////////////////////////////////////////////////////////////////////////////
void compute_sound_speed( const auto & eos,
                          const auto & density, 
                          const auto & internal_energy,
                          auto & sound_speed )
{
  auto tup = utils::zip( density, internal_energy, sound_speed );
  for ( auto it : tup ) {
    auto & d = it.template get<0>();
    auto & e = it.template get<1>();
    auto & s = it.template get<2>();
    s = eos.compute_sound_speed( d, e );
  }
}



////////////////////////////////////////////////////////////////////////////////
//! \brief comput the temperature
//!
//! This version is more efficient since there is only one virtual call to
//! evaluate a list of temperatures.
//! 
//! \param[in] density the density
//! \param[in] internal_energy the internal energy
//! \param[out] temperature the temperature
////////////////////////////////////////////////////////////////////////////////
void compute_temperature( const auto & eos,
                          const auto & density, 
                          const auto & internal_energy,
                          auto & temperature )
{
  auto tup = utils::zip( density, internal_energy, temperature );
  for ( auto it : tup ) {
    auto & d = it.template get<0>();
    auto & e = it.template get<1>();
    auto & t = it.template get<2>();
    t = eos.compute_temperature( d, e );
  }
}


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
