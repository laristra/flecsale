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

#include <flecsale/eos/ideal_gas.h>
#include <flecsale/mesh/burton/burton.h>
#include <flecsale/utils/lua_utils.h>

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

#ifdef HAVE_LUA

  //===========================================================================
  //! \brief Load the lua input file
  //! \param [in] file  The name of the lua file to load.
  //===========================================================================
  static auto load_lua(const std::string & file) 
  {
    // setup the python interpreter
    auto lua_state = flecsale::utils::lua_t();
    // load the test file
    lua_state.loadfile( file );

    // return the state
    return lua_state;
  }

#endif // HAVE_LUA

};

} // namespace hydro
} // namespace apps
