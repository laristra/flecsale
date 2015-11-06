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
 * \file eos_base.h
 * 
 * \brief Base class for the equation of state module.
 *
 ******************************************************************************/
#pragma once

//! system includes
#include <vector>

//! user includes
#include "ale/common/types.h"

namespace ale {
namespace eos {

////////////////////////////////////////////////////////////////////////////////
//! \brief Base EOS class
////////////////////////////////////////////////////////////////////////////////
class eos_base {


public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the real precision type
  using real_t = common::real_t;

#if 0


  //! \brief the index type
  using index_t = common::index_t;


  //! \brief define a function for accessing collections of data
  using container_t = std::vector<real_t>;

  //============================================================================
  // Constructors / Destructors
  //============================================================================


  //! \brief default constructor
  eos_base() {};

  //! \brief the destructor
  virtual ~eos_base() {};

  //! \brief disallow copy constructor
  eos_base(const eos_base&) = delete;
  
  //! \brief Disallow copying
  eos_base& operator=(const eos_base&) = delete;


  //============================================================================
  // Public Member Functions
  //============================================================================

  //! \brief Main driver to setup the eos using input data
  //! \param[in] input the input object with the data
  //virtual void initialize( const common::NodeType &input ) = 0;

  //! \brief Main driver to setup the eos object
  virtual void setup( void ) = 0;


  //! \brief return the density
  //! \return the density
  virtual real_t get_density( void ) const = 0;

  //! \brief return the internal energy
  //! \return the internal energy
  virtual real_t get_internal_energy( void ) const = 0;


  //! \brief compute the internal energy
  //!
  //! This version is less efficient since there is one virtual call for
  //! every internal energy evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] pressure the pressure
  //! \return the internal energy
  virtual real_t compute_internal_energy( real_t density, 
                                          real_t pressure ) const = 0;

  //! \brief compute the internal energy
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of internal energies.
  //! 
  //! \param[in]  density the density
  //! \param[in]  pressure the pressure
  //! \param[out] internal_energy the internal energy
  virtual void compute_internal_energy( const container_t & density, 
                                        const container_t & pressure,
                                        container_t & internal_energy ) const = 0;

  //! \brief compute the pressure
  //!
  //! This version is less efficient since there is one virtual call for
  //! every internal pressure evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the pressure
  virtual real_t compute_pressure( real_t density, 
                                   real_t internal_energy ) const = 0;

  //! \brief compute the pressure
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of pressures.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \param[out] pressure the pressure
  virtual void compute_pressure( const container_t & density, 
                                 const container_t & internal_energy,
                                 container_t & pressure ) const = 0;


  //! \brief comput the sound speed
  //!
  //! This version is less efficient since there is one virtual call for
  //! every sound speed evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the sound speed
  virtual real_t compute_sound_speed( real_t density, 
                                      real_t internal_energy ) const = 0;

  //! \brief comput the sound speed
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of sound speeds.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \param[out] sound_speed the sound speed
  virtual void compute_sound_speed( const container_t & density, 
                                    const container_t & internal_energy,
                                    container_t & sound_speed ) const = 0;

  //! \brief comput the temperature
  //!
  //! This version is less efficient since there is one virtual call for
  //! every temperature evaluation.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the temperature
  virtual real_t compute_temperature( real_t density, 
                                      real_t internal_energy ) const = 0;

  //! \brief comput the temperature
  //!
  //! This version is more efficient since there is only one virtual call to
  //! evaluate a list of temperatures.
  //! 
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \param[out] temperature the temperature
  virtual void compute_temperature( const container_t & density, 
                                    const container_t & internal_energy,
                                    container_t & temperature ) const = 0;

#endif

};

} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
