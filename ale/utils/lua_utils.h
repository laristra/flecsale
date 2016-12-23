/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utilities for using lua.
////////////////////////////////////////////////////////////////////////////////

#if HAVE_LUA 

// user includes
#include "ale/utils/string_utils.h"
#include "ale/utils/errors.h"

// use lua
extern "C" {
  #include <lua.h>
  #include <lualib.h>
  #include <lauxlib.h>
}

// system libraries
#include <sstream>

namespace ale {
namespace utils {

namespace detail {

template < typename T >
struct result {};

template <>
struct result<int> 
{
  static int get(lua_State * s) { return lua_tointeger(s,-1); }
};

template <>
struct result<double> 
{
  static int get(lua_State * s) { return lua_tonumber(s,-1); }
};

}

class lua {

  lua_State * state_ = nullptr;

public:

  lua(bool with_system = true) : state_(luaL_newstate()) 
  {
    if ( !state_ )
      raise_runtime_error("Cannot initialize lua state.");
    // open all system libraries
    if ( with_system )
      luaL_openlibs(state_);
  }

  ~lua() { lua_close(state_); }

  std::string get_stack(int i) const
  {
    std::stringstream os;
    auto t = lua_type(state_, i);
    switch (t) {
      case LUA_TSTRING:  // strings
        os << lua_tostring(state_, i);;
      case LUA_TBOOLEAN:  // booleans
        os << lua_toboolean(state_, i) ? "true" : "false";
      case LUA_TNUMBER:  // numbers
        os << lua_tonumber(state_, i);
      default:  // other values
        os << lua_typename(state_, t);
    }
    return os.str();
  }

  friend std::ostream& operator<<(std::ostream& os, const lua & s)  
  {  
    auto top = lua_gettop(s.state_);
    for (int i = 1; i <= top; i++) {  /* repeat for each level */
      os << s.get_stack(i);
      os << "  ";  // put a separator
    }
    os << std::endl; // end the listing
    return os;
  }

  void print_last_error()
  {
    std::cerr << get_stack(-1) << std::endl;
    lua_pop(state_, 1);
  }

  void run_string( std::string script )
  {
    auto ret = luaL_dostring(state_,script.c_str());
    if ( ret ) {
      print_last_error();
      raise_runtime_error("Cannot load buffer.");
    }
  }

  void loadfile( std::string file )
  {
    auto ret = luaL_dofile(state_,file.c_str());
    if ( ret ) {
      print_last_error();
      raise_runtime_error("Cannot load file.");
    }
  }

template<
  typename T,
  typename = std::enable_if_t< std::is_arithmetic<T>::value >
>
void push_arg(T && arg)
{
  lua_pushnumber(state_, std::forward<T>(arg) );
}


template< typename Arg >
void push_args( Arg&& arg )
{
  // push the last argument
  push_arg( std::forward<Arg>(arg) );
}

template< typename Arg, typename... Args >
void push_args( Arg&& arg, Args&&... args ) {
  // chop off an argument and push it
  push_arg( std::forward<Arg>(arg) );
  // set the remaining arguments
  push_args(std::forward<Args>(args)...);
}

  template <typename T, typename...Args> 
  T call_function(std::string name, Args&&...args)
  {
    // the function name
    lua_getglobal(state_, name.c_str());
    if ( !lua_isfunction(state_, -1) ) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" is not a function.");
    }
    // add the arguments
    push_args(std::forward<Args>(args)...);
    // call the function
    auto ret = lua_pcall(state_, sizeof...(args), 1, 0);
    if (ret) {
      print_last_error();
      raise_runtime_error("Problem calling \"" << name << "\".");
    }
    // make sure the result is non nill
    if (lua_isnil(state_, -1)) {
      print_last_error();
      raise_runtime_error("No result from function \"" << name << "\".");
    }
    // get the result
    auto res = detail::result<T>::get(state_);
    lua_pop(state_, 1);
    // return the result
    return res;
  }

};

#endif
  
} // namespace utils
} // namespace ale