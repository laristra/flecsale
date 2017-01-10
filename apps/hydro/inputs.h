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
#include <ale/mesh/burton/burton.h>
#include <ale/mesh/factory.h>
#include <ale/utils/lua_utils.h>

// system includes
#include <string>

namespace apps {
namespace hydro {

///////////////////////////////////////////////////////////////////////////////
//! \brief A struct that contains all the inputs.
///////////////////////////////////////////////////////////////////////////////
class inputs_t {
public:

  //! The number of dimensions
  static constexpr auto num_dimensions = 3;
 
  //! the mesh type
  using mesh_t = ale::mesh::burton_mesh_t<num_dimensions>; 
  //! the size type
  using size_t = typename mesh_t::size_t;
  //! the real type
  using real_t = typename mesh_t::real_t;
  //! the vector type
  using vector_t = typename mesh_t::vector_t;

  //! the ics function type
  //! \{
  using ics_return_t = std::tuple<real_t,vector_t,real_t>;
  using ics_function_t = std::function< ics_return_t(const vector_t & x) >;
  //! \}

  //! the mesh function type
  using mesh_function_t = std::function< mesh_t(void) >;

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

  //! \brief this is a lambda function to set the initial conditions
  static ics_function_t ics;

  //! \brief This function builds and returns a mesh
  static mesh_function_t make_mesh  ; 

  //! \brief Load the input file
  static void load(const std::string & file) 
  {
    auto ext = ale::utils::file_extension(file);
    if ( ext == "lua" ) {
      load_lua(file);
    }
    else
      raise_runtime_error(
        "Unknown file extension for \""<<file<<"\""
      );
  }
  
  //===========================================================================
  //! \brief Load the lua input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
  static void load_lua(const std::string & file) 
  {
    // setup the python interpreter
    auto state = ale::utils::lua_t();
    // load the test file
    state.loadfile( file );

    // get the hydro table
    auto hydro_ics = state["hydro"];

    // now set some inputs
    prefix = hydro_ics["prefix"].as<std::string>();
    postfix = hydro_ics["postfix"].as<std::string>();
    output_freq = hydro_ics["output_freq"].as<size_t>();
    CFL = hydro_ics["CFL"].as<real_t>();
    final_time = hydro_ics["final_time"].as<real_t>();
    max_steps = hydro_ics["max_steps"].as<size_t>();
    auto ics_func = hydro_ics["ics"];

    // set the ics function
    ics = [ics_func]( const vector_t & x )
      {
        real_t d, p;
        vector_t v(0);
        std::tie(d, v, p) = 
          ics_func(x[0], x[1], x[2]).as<real_t, vector_t, real_t>();
        return std::make_tuple( d, v, p );
      };
      
    // now set the mesh building function
    auto mesh_input = hydro_ics["mesh"];
    auto mesh_type = lua_try_access(mesh_input, "type", std::string );
    if ( mesh_type == "box" ) {
      auto dims = lua_try_access( mesh_input, "dimensions", std::vector<int> );
      auto lens = lua_try_access( mesh_input, "lengths", std::vector<real_t> );
      make_mesh = [=](void)
      {
        return ale::mesh::box<mesh_t>( 
          dims[0], dims[1], dims[2], lens[0], lens[1], lens[2] 
        );
      };
    }
    else if (mesh_type == "read" ) {
      auto file = lua_try_access( mesh_input, "file", std::string );
      make_mesh = [=](void)
      {
        mesh_t m;
        ale::mesh::read_mesh(file, m);
        return m;
      };
    }
    else {
      raise_implemented_error("Unknown mesh type \""<<mesh_type<<"\"");
    }
  }
};

} // namespace hydro
} // namespace apps
