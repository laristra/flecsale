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
#include "../inputs.h"

#include <flecsale/mesh/burton/burton.h>
#include <flecsale/mesh/factory.h>
#include <flecsale/utils/lua_utils.h>

// system includes
#include <string>

namespace apps {
namespace hydro {

///////////////////////////////////////////////////////////////////////////////
//! \brief A struct that contains all the inputs for a 3d case.
///////////////////////////////////////////////////////////////////////////////
class inputs_t : public inputs_<3> {
public:

  //! \brief The base class.
  using base_t = inputs_<3>;

  //===========================================================================
  //! \brief Load the input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
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
#ifdef HAVE_LUA
    // setup the lua interpreter and read the common inputs
    auto lua_state = base_t::load_lua(file);

    // get the hydro table
    auto hydro_input = lua_try_access( lua_state, "hydro" );

    // now set some dimension specific inputs

    // set the ics function
    auto ics_func = lua_try_access( hydro_input, "ics" );
    ics = 
      [ics_func]( const vector_t & x, const real_t & t )
      {
        real_t d, p;
        vector_t v(0);
        std::tie(d, v, p) = 
          ics_func(x[0], x[1], x[2], t).as<real_t, vector_t, real_t>();
        return std::make_tuple( d, std::move(v), p );
      };
      
    // now set the mesh building function
    auto mesh_input = hydro_input["mesh"];
    auto mesh_type = lua_try_access_as(mesh_input, "type", std::string );

    if ( mesh_type == "box" ) {
      auto dims = lua_try_access_as( mesh_input, "dimensions", array_t<int> );
      auto xmin = lua_try_access_as( mesh_input, "xmin", array_t<real_t> );
      auto xmax = lua_try_access_as( mesh_input, "xmax", array_t<real_t> );
      make_mesh = [dims,xmin,xmax](const real_t &)
      {
        return ale::mesh::box<mesh_t>( 
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
        ale::mesh::read_mesh(file, m);
        return m;
      };
    }
    else {
      raise_implemented_error("Unknown mesh type \""<<mesh_type<<"\"");
    }

    // now clear and reset the boundary conditions
    auto bcs_input = lua_try_access( hydro_input, "bcs" );
    bcs.clear();
    
    for ( int i=0; i<bcs_input.size(); ++i ) {
      // get each bc pair
      auto bc_input = bcs_input[i+1];
      auto bc_type = lua_try_access_as( bc_input, "type", std::string );
      auto bc_func = lua_try_access( bc_input, "func" );
      // make the boundary condition function
      auto bc_predicate = [=]( const vector_t & x, const real_t & t )
        { 
          return bc_func(x[0], x[1], x[2], t).as<bool>();
        };
      // make a new boundary condition type
      auto bc_object = bcs_ptr_t( 
        make_boundary_condition<num_dimensions>(bc_type) 
      );
      // now make the boundary condition pair
      auto new_bc = 
        std::make_pair( std::move(bc_object), std::move(bc_predicate) );
      bcs.emplace_back( std::move(new_bc) );
    }
    
#else

    raise_implemented_error( 
      "You need to link with lua in order to use lua functionality."
    );

#endif // HAVE_LUA
  }
};

} // namespace hydro
} // namespace apps
