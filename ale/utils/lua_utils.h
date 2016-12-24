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
#include <memory>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ale {
namespace utils {

using lua_state_ptr_t = std::shared_ptr<lua_State>;
 
template < typename T >
struct lua_value {};

template <>
struct lua_value<long long> 
{
  static long long get(lua_State * s, int index) 
  {
    auto t = lua_type(s, index);
    if ( t != LUA_TNUMBER )
     raise_runtime_error( "Invalid conversion of type \"" <<
      lua_typename(s, t) << "\" to int."
    );
    return lua_tointeger(s, index);
  }
};

template <>
struct lua_value<unsigned long long> 
{
  static unsigned long long get(lua_State * s, int index) 
  {
    return static_cast<unsigned long long>(lua_value<long long>::get(s, index));
  }
};

template <>
struct lua_value<int> 
{
  static int get(lua_State * s, int index) 
  {
    return static_cast<int>(lua_value<long long>::get(s, index));
  }
};

template <>
struct lua_value<unsigned int> 
{
  static unsigned int get(lua_State * s, int index) 
  {
    return static_cast<unsigned int>(lua_value<long long>::get(s, index));
  }
};

template <>
struct lua_value<double> 
{
  static double get(lua_State * s, int index) 
  {
    auto t = lua_type(s, index);
    if ( t != LUA_TNUMBER )
     raise_runtime_error( "Invalid conversion of type \"" <<
      lua_typename(s, t) << "\" to double."
    );
    return lua_tonumber(s, index);
  }
};

template <>
struct lua_value<float> 
{
  static float get(lua_State * s, int index) 
  {
    return static_cast<float>(lua_value<double>::get(s, index));
  }
};

template <>
struct lua_value<bool> 
{
  static bool get(lua_State * s, int index) 
  {
    return lua_toboolean(s, index);
  }
};

template <>
struct lua_value<std::string> 
{
  static std::string get(lua_State * s, int index) 
  {
    return lua_tostring(s, index);
  }
};

class lua_base_t {

protected:

  lua_state_ptr_t state_;

public:

  lua_base_t()
    : state_( luaL_newstate(), [](lua_State * s) { lua_close(s); } )
  {}

  lua_base_t(const lua_state_ptr_t & state) 
    : state_(state)
  {}
  
  auto state() { return state_.get(); }
  auto state() const { return state_.get(); }

  std::string get_row(int i) const
  {
    auto s = state();
    std::stringstream os;
    auto t = lua_type(s, i);
    switch (t) {
      case LUA_TSTRING:  // strings
        os << lua_tostring(s, i);;
      case LUA_TBOOLEAN:  // booleans
        os << lua_toboolean(s, i) ? "true" : "false";
      case LUA_TNUMBER:  // numbers
        os << lua_tonumber(s, i);
      default:  // other values
        os << lua_typename(s, t);
    }
    return os.str();
  }
  
  std::ostream& dump_stack(std::ostream& os) const 
  {
    auto s = state();
    auto top = lua_gettop(s);
    for (int i = 1; i <= top; i++) {  /* repeat for each level */
      os << get_row(i);
      os << "  ";  // put a separator
    }
    os << std::endl; // end the listing
    return os;
  }
  
  void print_last_error()
  {
    std::cerr << get_row(-1) << std::endl;
    lua_pop(state(), 1);
  }

  friend std::ostream& operator<<(std::ostream& os, const lua_base_t & s)  
  { 
    return s.dump_stack(os); 
  }
};

class lua_result_t : public lua_base_t {

  int index_;
  int num_args_;

  template< typename Tup, int I, typename Arg >
  void unpack_args( Tup & tup ) const 
  {
    std::get<I>(tup) = lua_value<Arg>::get( state(), index_-num_args_+I+1 );
    // stop recursion
  }
  
  template< 
    typename Tup, int I, typename Arg1, typename Arg2, typename... Args
  >
  void unpack_args( Tup & tup ) const 
  {
    std::get<I>(tup) = lua_value<Arg1>::get( state(), index_-num_args_+I+1 );
    // recursively extract tuple 
    unpack_args<Tup,I+1,Arg2,Args...>( tup );
  }

  void check_args( int n ) const
  {
    if ( n > num_args_ ) {
      raise_runtime_error(
        "Requested \""<<n<<"\" arguments, but function " <<
        "return \"" << num_args_ << "\"."
      );
    }
  }

public:

  lua_result_t() = delete;

  lua_result_t(lua_state_ptr_t & state, int index, int num_args) 
    : lua_base_t(state), index_(index), num_args_(num_args)
  {}

  template< typename T >
  explicit operator T() const {
    check_args(1);
    return lua_value<T>::get(state(), index_);
  }

  template< typename...Args >
  explicit operator std::tuple<Args...>() const {
    constexpr int N = sizeof...(Args);
    check_args(N);
    using Tup = std::tuple<Args...>; 
    Tup tup;
    unpack_args<Tup,0,Args...>(tup);
    return tup;
  }

  template< typename T >
  auto as() const {
    return static_cast<T>(*this);
  }
  
  template< typename T1, typename T2, typename...Ts >
  auto as() const {
    return static_cast< std::tuple<T1,T2,Ts...> >(*this);
  }
};

class lua_global_t : public lua_base_t {

  std::string name_;
  int index_;

  template< typename T >
  void push_arg(T && arg)
  {
    auto s = state();
    if ( std::is_integral<T>::value )
      lua_pushinteger(s, std::forward<T>(arg) );
    else if ( std::is_arithmetic<T>::value )
      lua_pushnumber(s, std::forward<T>(arg) );
    else
      raise_runtime_error(
        "Trying to push unknown type \"" << typeid(T).name() <<
        "\" onto stack."
      );
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

public:

  lua_global_t() = delete;

  lua_global_t(lua_state_ptr_t & state, const std::string & name, int index) 
    : lua_base_t(state), name_(name), index_(index)
  {}
  
  template < typename...Args >
  auto operator()(Args&&...args) 
  {
    auto s = state();
    // get the current stack position ( the function and arguments will get
    // deleted from the stack after the function call )
    auto pos0 = lua_gettop(s);
    // push the function onto the stack
    lua_pushvalue( s, index_ );
    // now make sure its a function
    if ( !lua_isfunction(s, -1) ) {
      print_last_error();
      raise_runtime_error("\"" << name_ << "\" is not a function.");
    }
    // add the arguments
    push_args(std::forward<Args>(args)...);
    // call the function
    auto ret = lua_pcall(s, sizeof...(args), LUA_MULTRET, 0);
    if (ret) {
      print_last_error();
      raise_runtime_error("Problem calling \"" << name_ << "\".");
    }
    // make sure the result is non nill
    if (lua_isnil(s, -1)) {
      print_last_error();
      raise_runtime_error("No result from function \"" << name_ << "\".");
    }
    // get the result
    auto pos1 = lua_gettop(s);
    return lua_result_t( state_, pos1, pos1-pos0 );
   }
};


class lua_t : public lua_base_t {

public:

  lua_t(bool with_system = true) : lua_base_t() 
  {
    if ( !state_ )
      raise_runtime_error("Cannot initialize lua state.");
    // open all system libraries
    if ( with_system )
      luaL_openlibs(state());
  }

  void run_string( const std::string & script )
  {
    auto ret = luaL_dostring(state(),script.c_str());
    if ( ret ) {
      print_last_error();
      raise_runtime_error("Cannot load buffer.");
    }
  }

  void loadfile( const std::string & file )
  {
    auto ret = luaL_dofile(state(),file.c_str());
    if ( ret ) {
      print_last_error();
      raise_runtime_error("Cannot load file.");
    }
  }

  template <typename...Args> 
  auto operator[]( const std::string & name )
  {
    auto s = state();
    // the function name
    lua_getglobal(s, name.c_str());
    // return the global object with a pointer to a location in the stack
    return lua_global_t( state_, name, lua_gettop(s) );
  }

};

#endif
  
} // namespace utils
} // namespace ale
