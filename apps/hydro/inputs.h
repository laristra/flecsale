/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines a struct that contains all the inputs.
///////////////////////////////////////////////////////////////////////////////
#pragma once 

// user includes
#include <flecsale/eos/eos_base.h>
#include <flecsale/eos/ideal_gas.h>
#include <flecsale/mesh/burton/burton.h>
#include <flecsale/utils/lua_utils.h>

// system includes
#include <string>

namespace apps {
namespace hydro {

///////////////////////////////////////////////////////////////////////////////
//! \brief A struct that contains all the inputs for a the base case.
///////////////////////////////////////////////////////////////////////////////
template< int N >
class inputs_ {
public:

  //! The number of dimensions
  static constexpr auto num_dimensions = N;
 
  //! the mesh type
  using mesh_t = ale::mesh::burton_mesh_t<num_dimensions>; 
  //! the size type
  using size_t = typename mesh_t::size_t;
  //! the real type
  using real_t = typename mesh_t::real_t;
  //! the vector type
  using vector_t = typename mesh_t::vector_t;
  //! the eos type
  using eos_t = ale::eos::eos_base_t<real_t>;

  //! a dimensioned array type helper
  template< typename T>
  using array_t = std::array<T, num_dimensions>;

  //! the ics function type
  //! \{
  using ics_return_t = std::tuple<real_t,vector_t,real_t>;
  using ics_function_t = 
    std::function< ics_return_t(const vector_t & x, const real_t & t) >;
  //! \}

  //! the mesh function type
  using mesh_function_t = std::function< mesh_t(const real_t & t) >;

  //! \brief the case prefix and postfix
  //! \{
  static std::string prefix;
  static std::string postfix;
  //! \}

  //! \brief output frequency
  static size_t output_freq;

  //! \brief the CFL and final solution time
  //! \{
  static real_t CFL;
  static real_t final_time;
  static size_t max_steps;
  //! \}

  //! \brief the equation of state
  static std::shared_ptr<eos_t> eos;

  //! \brief this is a lambda function to set the initial conditions
  static ics_function_t ics;

  //! \brief This function builds and returns a mesh
  static mesh_function_t make_mesh; 

#ifdef HAVE_LUA

  //===========================================================================
  //! \brief Load the lua input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
  static auto load_lua(const std::string & file) 
  {
    // setup the python interpreter
    auto lua_state = ale::utils::lua_t();
    // load the test file
    lua_state.loadfile( file );

    // get the hydro table
    auto hydro_input = lua_try_access( lua_state, "hydro" );

    // now set some inputs
    prefix = lua_try_access_as( hydro_input, "prefix", std::string );
    postfix = lua_try_access_as( hydro_input, "postfix", std::string );
    output_freq = lua_try_access_as( hydro_input, "output_freq", size_t );
    CFL = lua_try_access_as( hydro_input, "CFL", real_t );
    final_time = lua_try_access_as( hydro_input, "final_time", real_t );
    max_steps = lua_try_access_as( hydro_input, "max_steps", size_t );

    // setup the equation of state
    auto eos_input = lua_try_access( hydro_input, "eos" );
    auto eos_type = lua_try_access_as( eos_input, "type", std::string );
    if ( eos_type == "ideal_gas" ){
      using ideal_gas_t = ale::eos::ideal_gas_t<real_t>;
      auto g  = lua_try_access_as( eos_input, "gas_constant", real_t );
      auto cv = lua_try_access_as( eos_input, "specific_heat", real_t );
      eos = std::make_shared<ideal_gas_t>( g, cv );
    }
    else {
      raise_implemented_error("Unknown eos type \""<<eos_type<<"\"");
    }

    // return the state
    return lua_state;
  }

#endif // HAVE_LUA

};

} // namespace hydro
} // namespace apps
