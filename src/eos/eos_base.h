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
 * \file Eos.h
 * 
 * \brief Base class for the equation of state module.
 *
 ******************************************************************************/
#pragma once

//! system includes

//! user includes
#include "common/types.h"

namespace ale {
namespace eos {

////////////////////////////////////////////////////////////////////////////////
//! \brief Base EOS class
////////////////////////////////////////////////////////////////////////////////
class eos_base {


public:

  //============================================================================
  // Constructors / Destructors
  //============================================================================

  //! \brief the real precision type
  using real_t = common::real_t;

  //============================================================================
  // Constructors / Destructors
  //============================================================================


  //! \brief default constructor
  eos_base() {};

  //! \brief the destructor
  virtual ~eos_base() {};

  //! \brief disallow copy constructor
  Eos(const Eos&) = delete;
  
  //! \brief Disallow copying
  Eos& operator=(const Eos&) = delete;


  //============================================================================
  // Public Member Functions
  //============================================================================

  //! \brief Main driver to setup the eos using input data
  //! \param[in] input the input object with the data
  virtual void initialize( const common::NodeType &input ) = 0;

  //! \brief Main driver to setup the eos object
  virtual void setup( void ) = 0;


  //! \brief return the density
  //! \return the density
  virtual real_t get_density( void ) const = 0;

  //! \brief return the internal energy
  //! \return the internal energy
  virtual real_t get_internal_energy( void ) const = 0;


  //! \brief compute the internal energy
  //! \param[in] density the density
  //! \param[in] pressure the pressure
  //! \return the internal energy
  virtual real_t compute_internal_energy( real_t density, 
                                          real_t pressure ) const = 0;

  //! \brief compute the pressure
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the pressure
  virtual real_t computePressure( real_t density, 
                                  real_t internal_energy ) const = 0;


  //! \brief comput the sound speed
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the sound speed
  virtual real_t computeSoundSpeed( real_t density, 
                                    real_t internal_energy ) const = 0;


  //! \brief comput the temperature
  //! \param[in] density the density
  //! \param[in] internal_energy the internal energy
  //! \return the temperature
  virtual real_t computeTemperature( real_t density, 
                                     real_t internal_energy ) const = 0;


};

} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
