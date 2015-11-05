/********************************************************************
      ___           ___           ___           ___       ___     
     /\  \         /\  \         /\  \         /\__\     /\  \    
    /::\  \       /::\  \       /::\  \       /:/  /    /::\  \   
   /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/  /    /:/\:\  \  
  /::\~\:\  \   /::\~\:\  \   /::\~\:\  \   /:/  /    /::\~\:\  \ 
 /:/\:\ \:\__\ /:/\:\ \:\__\ /:/\:\ \:\__\ /:/__/    /:/\:\ \:\__\
 \/_|::\/:/  / \:\~\:\ \/__/ \/__\:\/:/  / \:\  \    \:\~\:\ \/__/
    |:|::/  /   \:\ \:\__\        \::/  /   \:\  \    \:\ \:\__\  
    |:|\/__/     \:\ \/__/        /:/  /     \:\  \    \:\ \/__/  
    |:|  |        \:\__\         /:/  /       \:\__\    \:\__\    
     \|__|         \/__/         \/__/         \/__/     \/__/    

********************************************************************/
/*!
 *
 * \file IdealGas.cpp
 * 
 * \brief Ideal gas EOS class implementation.
 *
 ********************************************************************/

//! main include
#include "physics/IdealGas.h"

//! system includes
#include <cassert>
#include <iostream>
#include <string>

//! user includes
#include "errors/ErrorUtils.h"

namespace ReALE {
namespace physics {


  ///////////////////////////////////////////////////////////////////
  // The default constructor
  ///////////////////////////////////////////////////////////////////
  IdealGas::IdealGas() : 
    Gamma_(1.4), Density_(1.0), SpecificHeatV_(1.0), InternalEnergy_(1.0)
  {
  };

  ///////////////////////////////////////////////////////////////////
  // The destructor
  ///////////////////////////////////////////////////////////////////
  IdealGas::~IdealGas (void) 
  { };
  
  ///////////////////////////////////////////////////////////////////
  // Main driver to setup the eos using input data
  ///////////////////////////////////////////////////////////////////
  void IdealGas::initialize( const common::NodeType &input ) 
  {  
    if ( input["gamma"] ) {
      Gamma_ = input["gamma"].as<common::RealType>();
      if ( Gamma_ < 0.0 )
        raise_logic_error("Negative gamma used in material input");
    }

    if ( input["density"] ) {
      Density_ = input["density"].as<common::RealType>();  
      if ( Density_ < 0.0 )
        raise_logic_error("Negative density used in material input");
    }

    if ( input["specific_heat"] ) {
      SpecificHeatV_  = input["specific_heat"].as<common::RealType>();
      if ( SpecificHeatV_ < 0.0 )
        raise_logic_error("Negative specific heat used in material input");
    }

    // initialize state with temperature
    if ( input["temperature"] ) {
      common::RealType Temperature = input["temperature"].as<common::RealType>();  
      if ( Temperature < 0.0 )
        raise_logic_error("Negative temperature used in material input");      
      InternalEnergy_ = SpecificHeatV_ * Temperature;
    }

    // initialize state with pressure
    if ( input["pressure"] ) {
      common::RealType Pressure = input["pressure"].as<common::RealType>();  
      if ( Pressure < 0.0 )
        raise_logic_error("Negative temperature used in material input");      
      InternalEnergy_ = Pressure / ( Density_*(Gamma_-1.0) );
    }

    if ( input["temperature"] && input["pressure"] ) {
      std::string name = input["name"].as<std::string>("IdealGas");
      std::cout << std::endl
                << "Warning: Both a temperature and pressure were specified for material <" 
                << name << ">" 
                << std::endl
                << std::endl;
    }

  }

  ///////////////////////////////////////////////////////////////////
  // setup the ideal gas
  ///////////////////////////////////////////////////////////////////
  void IdealGas::setup(void)
  { 
  }


  ///////////////////////////////////////////////////////////////////
  // get the internal energy
  ///////////////////////////////////////////////////////////////////
  common::RealType IdealGas::getInternalEnergy( void ) const
  { return InternalEnergy_; }


  ///////////////////////////////////////////////////////////////////
  // get the density
  ///////////////////////////////////////////////////////////////////
  common::RealType IdealGas::getDensity( void ) const
  { return Density_; }


  ///////////////////////////////////////////////////////////////////
  // compute the internal energy
  ///////////////////////////////////////////////////////////////////
  common::RealType IdealGas::computeInternalEnergy( const common::RealType density, 
                                                    const common::RealType pressure ) const
  { return pressure / ( density*(Gamma_-1.0) ); }

  ///////////////////////////////////////////////////////////////////
  // compute the pressure
  ///////////////////////////////////////////////////////////////////
  common::RealType IdealGas::computePressure( const common::RealType density, 
                                              const common::RealType internal_energy ) const
  { return (Gamma_-1.0) * density * internal_energy; }


  ///////////////////////////////////////////////////////////////////
  // comput the sound speed
  ///////////////////////////////////////////////////////////////////
  common::RealType IdealGas::computeSoundSpeed( const common::RealType density, 
                                                const common::RealType internal_energy ) const
  { 
    assert( internal_energy>0.0 );
    return std::sqrt( Gamma_ * (Gamma_-1.0) * internal_energy );
  }



  ///////////////////////////////////////////////////////////////////
  // comput the temperature
  ///////////////////////////////////////////////////////////////////
  common::RealType IdealGas::computeTemperature( const common::RealType density, 
                                                 const common::RealType internal_energy ) const
  { return internal_energy / SpecificHeatV_; }


} // namespace
} // namespace
