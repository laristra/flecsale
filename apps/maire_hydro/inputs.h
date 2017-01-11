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
#include "types.h"

#include <ale/mesh/burton/burton.h>
#include <ale/utils/lua_utils.h>

// system includes
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
  using ics_function_t = 
    std::function< ics_return_t(const vector_t & x, const real_t & t) >;
  //! \}

  //! the mesh function type
  using mesh_function_t = std::function< mesh_t(const real_t & t) >;

  //! the bcs function type
  //! \{
  using bcs_t = boundary_condition_t<num_dimensions>;
  using bcs_ptr_t = std::shared_ptr< bcs_t >;
  using bcs_function_t = 
    std::function< bool(const vector_t & x, const real_t & t) >;
  using bcs_list_t = std::vector< std::pair< bcs_ptr_t, bcs_function_t > >;
  //! \}

  //! \brief the case prefix and postfix
  //! \{
  static std::string prefix;
  static std::string postfix;
  //! \}

  //! \brief output frequency
  static size_t output_freq;

  //! \brief the CFL and final solution time
  //! \{
  static time_constants_t CFL;
  static real_t final_time;
  static real_t initial_time_step;
  static size_t max_steps;
  //! \}

  //! \brief this is a lambda function to set the initial conditions
  static ics_function_t ics;

  //! \brief This function builds and returns a mesh
  static mesh_function_t make_mesh; 

  //! \brief this is a list of lambda functions to set the boundary conditions
  static bcs_list_t bcs;


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
    final_time = lua_try_access_as( hydro_input, "final_time", real_t );
    max_steps = lua_try_access_as( hydro_input, "max_steps", size_t );
    initial_time_step = lua_try_access_as( hydro_input, "initial_time_step", real_t );

    auto cfl_ics = lua_try_access( hydro_input, "CFL" );
    CFL.accoustic = lua_try_access_as( cfl_ics, "accoustic", real_t );
    CFL.volume    = lua_try_access_as( cfl_ics, "volume",    real_t );
    CFL.growth    = lua_try_access_as( cfl_ics, "growth",    real_t );

    // return the state
    return lua_state;
  }
};

} // namespace hydro
} // namespace apps
