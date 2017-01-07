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
  using ics_return_t = std::tuple<real_t,vector_t,real_t>;

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
  static ics_return_t ics(const vector_t & x);

  //! \brief This function builds and returns a mesh
  static mesh_t make_mesh(); 

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
  
  //! \brief Load the lua input file
  static void load_lua(const std::string & file) 
  {
    
  }
};

} // namespace hydro
} // namespace apps
