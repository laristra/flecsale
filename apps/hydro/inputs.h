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
#include <ale/utils/lua_utils.h>

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
  static mesh_function_t make_mesh; 

  //===========================================================================
  //! \brief Load the lua input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
  static auto load_lua(const std::string & file) 
  {
    // setup the python interpreter
    auto state = ale::utils::lua_t();
    // load the test file
    state.loadfile( file );

    // get the hydro table
    auto hydro_ics = state["hydro"];

    // now set some inputs
    prefix = lua_try_access( hydro_ics, "prefix", std::string );
    postfix = lua_try_access( hydro_ics, "postfix", std::string );
    output_freq = lua_try_access( hydro_ics, "output_freq", size_t );
    CFL = lua_try_access( hydro_ics, "CFL", real_t );
    final_time = lua_try_access( hydro_ics, "final_time", real_t );
    max_steps = lua_try_access( hydro_ics, "max_steps", size_t );

    // return the state
    return state;
  }
};

} // namespace hydro
} // namespace apps
