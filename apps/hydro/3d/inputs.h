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
#include <ristra/utils/string_utils.h>
#include "../types.h"

// system includes
#include <iomanip>
#include <string>

namespace apps {
namespace hydro {

///////////////////////////////////////////////////////////////////////////////
//! \brief A struct that contains all the inputs for a 3d case.
///////////////////////////////////////////////////////////////////////////////
class inputs_t {
public:
  
  //! the size type
  using size_t = mesh_t::size_t;
  //! the real type
  using real_t = mesh_t::real_t;
  //! the vector type
  using vector_t = mesh_t::vector_t;

  //! the number of dimensions
  static constexpr auto num_dimensions = mesh_t::num_dimensions;

  //! the eos type
  using eos_t = apps::hydro::eos_t;

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
  static eos_t eos;

  //! \brief this is a lambda function to set the initial conditions
  static ics_function_t ics;

  //===========================================================================
  //! \brief Load the input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
  static void load(const std::string & file)
  {
    auto ext = ristra::utils::file_extension(file);
    if ( ext == "lua" ) {
      load_lua(file);
    }
    else
      THROW_RUNTIME_ERROR(
        "Unknown file extension for \""<<file<<"\""
      );
  }
  
  //===========================================================================
  //! \brief Load the lua input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
  static void load_lua(const std::string & file) 
  {
#ifdef HAVE_LUA

    // setup the python interpreter
    auto lua_state = flecsale::utils::lua_t();
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
      using ideal_gas_t = flecsale::eos::ideal_gas_t<real_t>;
      auto g  = lua_try_access_as( eos_input, "gas_constant", real_t );
      auto cv = lua_try_access_as( eos_input, "specific_heat", real_t );
      eos = ideal_gas_t( g, cv );
    }
    else {
      THROW_IMPLEMENTED_ERROR("Unknown eos type \""<<eos_type<<"\"");
    }

    // now set some dimension specific inputs

    // set the ics function
    auto ics_func = lua_try_access( hydro_input, "ics" );
    ics = [ics_func]( const vector_t & x, const real_t & t )
      {
        real_t d, p;
        vector_t v(0);
        std::tie(d, v, p) = 
          ics_func(x[0], x[1], x[2], t).as<real_t, vector_t, real_t>();
        return std::make_tuple( d, std::move(v), p );
      };
      
    // now set the mesh building function
    auto mesh_input = lua_try_access( hydro_input, "mesh" );
    auto mesh_type = lua_try_access_as(mesh_input, "type", std::string );
    if ( mesh_type == "box" ) {
      auto dims = lua_try_access_as( mesh_input, "dimensions", array_t<int> );
      auto xmin = lua_try_access_as( mesh_input, "xmin", array_t<real_t> );
      auto xmax = lua_try_access_as( mesh_input, "xmax", array_t<real_t> );
      make_mesh = [dims,xmin,xmax](const real_t &)
      {
        return flecsale::mesh::box<mesh_t>( 
          dims[0], dims[1], dims[2], 
          xmin[0], xmin[1], xmin[2], 
          xmax[0], xmax[1], xmax[2]
        );
      };
    }
    else if (mesh_type == "read" ) {
      auto file = lua_try_access_as( mesh_input, "file", std::string );
      make_mesh = [file](const real_t &)
      {
        mesh_t m;
        //flecsale::mesh::read_mesh(file, m);
        return m;
      };
    }
    else {
      THROW_IMPLEMENTED_ERROR("Unknown mesh type \""<<mesh_type<<"\"");
    }

#else

    THROW_IMPLEMENTED_ERROR( 
      "You need to link with lua in order to use lua functionality."
    );

#endif // HAVE_LUA
  }
};

} // namespace hydro
} // namespace apps
