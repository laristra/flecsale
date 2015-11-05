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
 * \file IdealGas.h
 * 
 * \brief Ideal gas equation of state implementation.
 *
 ********************************************************************/
#ifndef IDEALGAS_H
#define IDEALGAS_H

//! Machine Config
#include "common/Config.h"

//! user includes
#include "physics/Eos.h"

namespace ReALE {
namespace physics {

  ///////////////////////////////////////////////////////////////////
  //! \brief Ideal Gas EOS class
  ///////////////////////////////////////////////////////////////////
  class IdealGas : public Eos {
    
  public:
    
    //===============================================================
    // Constructors / Destructors
    //===============================================================
    
    
    //! \brief default constructor
    IdealGas();
    
    //! \brief the destructor
    virtual ~IdealGas();
    
    
    //===============================================================
    // Public Member Functions
    //===============================================================
    
    
    //! \brief Main driver to setup the eos using input data
    //! \param[in] input the input object with the data
    virtual void initialize( const common::NodeType &input );

    //! \brief Main driver to setup the eos object
    virtual void setup( void );


    //! \brief return the density
    //! \return the density
    virtual common::RealType getDensity( void ) const;

    //! \brief return the internal energy
    //! \return the internal energy
    virtual common::RealType getInternalEnergy( void ) const;


    //! \brief compute the internal energy
    //! \param[in] density the density
    //! \param[in] pressure the pressure
    //! \return the internal energy
    virtual common::RealType computeInternalEnergy( const common::RealType density, 
                                                    const common::RealType pressure ) const;

    //! \brief compute the pressure
    //! \param[in] density the density
    //! \param[in] internal_energy the internal energy
    //! \return the pressure
    virtual common::RealType computePressure( const common::RealType density, 
                                              const common::RealType internal_energy ) const;


    //! \brief comput the sound speed
    //! \param[in] density the density
    //! \param[in] internal_energy the internal energy
    //! \return the sound speed
    virtual common::RealType computeSoundSpeed( const common::RealType density, 
                                                const common::RealType internal_energy ) const;


    //! \brief comput the temperature
    //! \param[in] density the density
    //! \param[in] internal_energy the internal energy
    //! \return the temperature
    virtual common::RealType computeTemperature( const common::RealType density, 
                                                 const common::RealType internal_energy ) const;
    
    
  protected:
    
    //===============================================================
    // Member variables
    //===============================================================
    
    //! the ideal gas specific heat ration
    common::RealType Gamma_; 

    //! the specific heat
    common::RealType SpecificHeatV_; 

    //! the reference density
    common::RealType Density_; 
    //! the reference internal energy
    common::RealType InternalEnergy_; 
    
private:
    
    //===============================================================
    // Deleted member functions
    //===============================================================
    
    //! \brief disallow copy constructor
    IdealGas(const IdealGas&);
    
    //! \brief Disallow copying
    IdealGas& operator=(const IdealGas&);


};

} // namespace
} // namespace


#endif // IDEALGAS_H


