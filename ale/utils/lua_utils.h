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
#include <iostream>
#include <iomanip>
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
        os << "string >> " << lua_tostring(s, i);
        break;
      case LUA_TBOOLEAN:  // booleans
        os << "bool   >> " << (lua_toboolean(s, i) ? "true" : "false");
        break;
      case LUA_TNUMBER:  // numbers
        os << "number >> " << lua_tonumber(s, i);
        break;
      default:  // other values
        os << "other  >> " << lua_typename(s, t);
        break;
    }
    return os.str();
  }
  
  std::ostream& dump_stack(std::ostream& os) const 
  {
    auto s = state();
    auto top = lua_gettop(s);
    os << "Row : Type   >> Value" << std::endl;
    for (int i = 1; i <= top; i++) {  /* repeat for each level */
      os << std::setw(3) << i << " : ";
      os << get_row(i);
      os << std::endl;  // put a separator
    }
    return os;
  }
  
  void print_last_error() const
  {
    auto s = state();
    auto top = lua_gettop(s);
    std::cerr << "Row : Type   >> Value" << std::endl;
    std::cerr << std::setw(3) << top << " : " << get_row(-1) << std::endl;
    //lua_pop(state(), 1);
  }

  friend std::ostream& operator<<(std::ostream& os, const lua_base_t & s)  
  { 
    return s.dump_stack(os); 
  }
};

class lua_result_t : public lua_base_t {

  std::string name_;
  int index_;

  template< typename Tup, int I, typename Arg >
  void unpack_args( Tup & tup ) const 
  {
    constexpr auto N = std::tuple_size<Tup>::value;
    std::get<I>(tup) = lua_value<Arg>::get( state(), index_-N+I+1 );
    // stop recursion
  }
  
  template< 
    typename Tup, int I, typename Arg1, typename Arg2, typename... Args
  >
  void unpack_args( Tup & tup ) const 
  {
    constexpr auto N = std::tuple_size<Tup>::value;
    std::get<I>(tup) = lua_value<Arg1>::get( state(), index_-N+I+1 );
    // recursively extract tuple 
    unpack_args<Tup,I+1,Arg2,Args...>( tup );
  }

  void push_arg(int i) const
  { lua_pushinteger( state(), i ); }
  
  void push_arg(long long i) const 
  { lua_pushinteger( state(), i ); }

  void push_arg(float x) const
  { lua_pushnumber( state(), x ); }

  void push_arg(double x) const
  { lua_pushnumber( state(), x ); }

  void push_arg(bool b) const
  { lua_pushboolean( state(), b ); }

  void push_arg(const char * s) const
  { lua_pushstring( state(), s ); }

  void push_arg(const std::string & s) const
  { lua_pushlstring( state(), s.c_str(), s.size() ); }
  
  void push_args() const
  {}
  
  template< typename Arg >
  void push_args( Arg&& arg ) const
  {
    // grow the stack
    check_stack(1);
    // push the last argument
    push_arg( std::forward<Arg>(arg) );
  }
  
  template< typename Arg, typename... Args >
  void push_args( Arg&& arg, Args&&... args ) const
  {
    // chop off an argument and push it
    push_arg( std::forward<Arg>(arg) );
    // set the remaining arguments
    push_args(std::forward<Args>(args)...);
  }
  
  void check_stack(int extra) const
  {
    auto s = state();
    auto ret = lua_checkstack(s, extra);
    if ( !ret ) {
      std::ostringstream ss;
      ss << "Cannot grow stack " << extra << " slots operating on element \""
         << name_ << "\"." << std::endl << "Current stack size is " 
         << lua_gettop(s) << ".";
      raise_runtime_error(ss.str());
    }
  } 

  void check_table(int index, const std::string & name) const
  {
    if ( !lua_istable(state(), index) ) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" is not a table.");
    }
  }

  void check_result(int index, const std::string & name) const
  {
    if (lua_isnil(state(), index)) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" returned nil.");
    }
  }

  void check_function(int index, const std::string & name) const
  {
    if ( !lua_isfunction(state(), index) ) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" is not a function.");
    }
  }

public:

  lua_result_t() = delete;

  lua_result_t(
    const lua_state_ptr_t & state, 
    const std::string & name, 
    int index
  ) : lua_base_t(state), name_(name), index_(index)
  {}

  template< typename T >
  explicit operator T() const {
    return lua_value<T>::get(state(), index_);
  }

  template< typename...Args >
  explicit operator std::tuple<Args...>() const {
    constexpr int N = sizeof...(Args);
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
  
  template < typename...Args >
  lua_result_t operator()(Args&&...args) const 
  {
    auto s = state();
    // get the current stack position ( the function and arguments will get
    // deleted from the stack after the function call )
    auto pos0 = lua_gettop(s);
    // push the function onto the stack
    check_stack(1);
    lua_pushvalue( s, index_ );
    // now make sure its a function
    check_function(index_, name_);
    // add the arguments
    push_args(std::forward<Args>(args)...);
    // keep track of the arguments
    auto pos1 = lua_gettop(s);
    std::stringstream args_ss;
    args_ss << "(";
    for ( auto p = pos0+1; p<=pos1; ++p ) {
      auto t = lua_type(s, p);
      args_ss << lua_typename(s, t);
      if (p<pos1) args_ss << ",";
    }
    args_ss << ")";
    // call the function
    auto ret = lua_pcall(s, sizeof...(args), LUA_MULTRET, 0);
    if (ret) {
      print_last_error();
      raise_runtime_error("Problem calling \"" << name_ << "\".");
    }
    // make sure the result is non nill
    check_result(-1,name_);
    // get the result
    pos1 = lua_gettop(s);
    return { state_, name_+args_ss.str(), pos1 };
  }

  lua_result_t operator[]( const std::string & name ) const &
  {
    auto s = state();
    auto new_name = name_+".[\""+name+"\"]";
    // make sure we are accessing a table
    check_table(index_, name_);
    // push the key onto the stack
    check_stack(1);
    lua_pushlstring(s, name.c_str(), name.size());
    // now get the table value
    lua_rawget(s, index_);
    // make sure returned result is not nil
    check_result(-1, new_name);
    // return the global object with a pointer to a location in the stack
    return { state_, new_name, lua_gettop(s) };
  }

  lua_result_t operator[]( int n ) const &
  {
    auto s = state();
    auto new_name = name_+"["+std::to_string(n)+"]";
    // make sure we are accessing a table
    check_table(index_, name_);
    // access the value
    lua_rawgeti(s, index_, n);
    // make sure returned result is not nil
    check_result(-1, new_name);
    // return the global object with a pointer to a location in the stack
    return { state_, new_name, lua_gettop(s) };
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

  bool run_string( const std::string & script )
  {
    auto ret = luaL_dostring(state(),script.c_str());
    if ( ret ) {
      print_last_error();
      raise_runtime_error("Cannot load buffer.");
    }
    return ret;
  }

  void loadfile( const std::string & file )
  {
    auto ret = luaL_dofile(state(),file.c_str());
    if ( ret ) {
      print_last_error();
      raise_runtime_error("Cannot load file.");
    }
  }

  lua_result_t operator[]( const std::string & name ) const &
  {
    auto s = state();
    // the function name
    lua_getglobal(s, name.c_str());
    // return the global object with a pointer to a location in the stack
    return { state_, name, lua_gettop(s) };
  }

  auto operator()( const std::string & script )
  {
    return run_string( script );
  }
};

#endif
  
} // namespace utils
} // namespace ale
