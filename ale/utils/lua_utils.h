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
#include <algorithm>
#include <cassert>
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <>
struct lua_value<long long> 
{
  static long long get(lua_State * s, int index = -1) 
  {
    if ( !lua_isnumber(s,index) )
     raise_runtime_error( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to int."
    );
    auto i = lua_tointeger(s, index);
    lua_remove(s,index);
    return i;
  }
};

template <>
struct lua_value<unsigned long long> 
{
  static unsigned long long get(lua_State * s, int index = -1) 
  {
    return static_cast<unsigned long long>(lua_value<long long>::get(s,index));
  }
};

template <>
struct lua_value<int> 
{
  static int get(lua_State * s, int index = -1) 
  {
    return static_cast<int>(lua_value<long long>::get(s,index));
  }
};

template <>
struct lua_value<unsigned int> 
{
  static unsigned int get(lua_State * s, int index = -1) 
  {
    return static_cast<unsigned int>(lua_value<long long>::get(s,index));
  }
};

template <>
struct lua_value<double> 
{
  static double get(lua_State * s, int index = -1) 
  {
    if ( !lua_isnumber(s,index) )
     raise_runtime_error( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to double."
    );
    auto x = lua_tonumber(s, index);
    lua_remove(s,index);
    return x;
  }
};

template <>
struct lua_value<float> 
{
  static float get(lua_State * s, int index = -1) 
  {
    return static_cast<float>(lua_value<double>::get(s, index));
  }
};

template <>
struct lua_value<bool> 
{
  static bool get(lua_State * s, int index = -1) 
  {
    if ( !lua_isboolean(s,index) )
     raise_runtime_error( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to bool."
    );
    auto b = lua_toboolean(s, index);
    lua_remove(s, index);
    return b;
  }
};

template <>
struct lua_value<std::string> 
{
  static std::string get(lua_State * s, int index = -1)
  {
    if ( !lua_isstring(s, index) )
     raise_runtime_error( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to string."
    );
    auto str = lua_tostring(s, index);
    lua_remove(s, index);
    return str;
  }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void lua_push(lua_State * s, int i)
{ lua_pushinteger( s, i ); }

void lua_push(lua_State * s, long long i)
{ lua_pushinteger( s, i ); }

void lua_push(lua_State * s, float x)
{ lua_pushnumber( s, x ); }

void lua_push(lua_State * s, double x) 
{ lua_pushnumber( s, x ); }

void lua_push(lua_State * s, bool b)
{ lua_pushboolean( s, b ); }

void lua_push(lua_State * s, const char * str)
{ lua_pushstring( s, str ); }

void lua_push(lua_State * s, const std::string & str) 
{ lua_pushlstring( s, str.c_str(), str.size() ); }
  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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
    if ( top ) {
      os << "Row : Type   >> Value" << std::endl;
      for (int i = 1; i <= top; i++) {  /* repeat for each level */
        os << std::setw(3) << i << " : ";
        os << get_row(i);
        os << std::endl;  // put a separator
      }
    }
    else {
      os << "(stack empty)" << std::endl;
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class lua_ref_t : lua_base_t {

  std::shared_ptr<int> ref_;

public:

  lua_ref_t() = delete;

  lua_ref_t ( const lua_state_ptr_t & state, int ref )
    : lua_base_t(state), 
      ref_( new int{ref},
        [s=state](int * r) 
        { luaL_unref(s.get(), LUA_REGISTRYINDEX, *r); } )
  {}

  lua_ref_t( const lua_state_ptr_t & state )
    : lua_ref_t(state, LUA_REFNIL)
  {}

  void push() const
  {
    lua_rawgeti(state(), LUA_REGISTRYINDEX, *ref_);
  }

};

lua_ref_t make_lua_ref(const lua_state_ptr_t & state)
{
  return { state, luaL_ref(state.get(), LUA_REGISTRYINDEX) };
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class lua_result_t : public lua_base_t {

  std::string name_;
  std::vector<lua_ref_t> refs_;

  template< typename Tup, int I >
  void get_results( Tup & tup ) const 
  {}
  
  template< 
    typename Tup, int I, typename Arg1, typename... Args
  >
  void get_results( Tup & tup ) const 
  {
    std::get<I>(tup) = lua_value<Arg1>::get( state() );
    // recursively extract tuple 
    get_results<Tup,I+1,Args...>( tup );
  }

  void push_args() const
  {}
  
  template< typename Arg, typename... Args >
  void push_args( Arg&& arg, Args&&... args ) const
  {
    // grow the stack
    check_stack(1);
    // chop off an argument and push it
    lua_push( state(), std::forward<Arg>(arg) );
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

  void check_table(const std::string & name) const
  {
    if ( !lua_istable(state(), -1) ) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" is not a table.");
    }
  }

  void check_result(const std::string & name) const
  {
    if (lua_isnil(state(), -1)) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" returned nil.");
    }
  }

  void check_function(const std::string & name) const
  {
    if ( !lua_isfunction(state(), -1) ) {
      print_last_error();
      raise_runtime_error("\"" << name << "\" is not a function.");
    }
  }

  void push_last() const
  {
    check_stack(1);
    refs_.back().push();
  }

  void push_all() const
  {
    // make sure the stack can handle this
    check_stack(refs_.size());
    // push everything onto the stack in reverse order
    for ( auto && r : refs_ ) std::forward<decltype(r)>(r).push();
  }

public:

  lua_result_t() = delete;

  lua_result_t(
    const lua_state_ptr_t & state, 
    const std::string & name, 
    lua_ref_t && ref
  ) : lua_base_t(state), name_(name)
  {
    refs_.emplace_back( std::move(ref) );
  }
  
  lua_result_t(
    const lua_state_ptr_t & state, 
    const std::string & name, 
    std::vector<lua_ref_t> && refs
  ) : lua_base_t(state), name_(name), refs_(std::move(refs))
  {}

  template< typename T >
  explicit operator T() const {
    if ( refs_.size() != 1 )
      raise_runtime_error( "Expecting 1 result, stack has " << refs_.size() );
    refs_.back().push();
    return lua_value<T>::get(state());
  }

  template< typename...Args >
  explicit operator std::tuple<Args...>() const {
    constexpr int N = sizeof...(Args);
    using Tup = std::tuple<Args...>; 
    if ( refs_.size() != N )
      raise_runtime_error( 
        "Expecting "<<N<<" results, stack has " << refs_.size() 
      );
    push_all();
    Tup tup;
    get_results<Tup,0,Args...>(tup);
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
    push_last();
    // now make sure its a function
    check_function(name_);
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
    check_result(name_);
    // figure out how much the stack grew
    pos1 = lua_gettop(s);
    auto num_results = pos1 - pos0;
    assert( num_results > 0 );
    // because there could be multiple results, get a reference to each one
    std::vector<lua_ref_t> ref_list;
    ref_list.reserve( num_results );
    for (int i=0; i<num_results; ++i)
      ref_list.emplace_back( std::move(make_lua_ref(state_)) );
    // get the result
    return { state_, name_+args_ss.str(), std::move(ref_list) };
  }

  lua_result_t operator[]( const std::string & name ) const &
  {
    auto s = state();
    auto new_name = name_+".[\""+name+"\"]";
    // push the table onto the stack
    push_last();
    // make sure we are accessing a table
    check_table(name_);
    // push the key onto the stack
    check_stack(1);
    lua_pushlstring(s, name.c_str(), name.size());
    // now get the table value, the key gets pushed from the stack
    lua_rawget(s, -2);
    // make sure returned result is not nil
    check_result(new_name);
    // get a reference to the result, and pop the table from the stack
    auto ref = make_lua_ref(state_);
    lua_pop(s, -1);
    // return the global object with a pointer to a location in the stack
    return { state_, new_name, std::move(ref) };
  }

  lua_result_t operator[]( int n ) const &
  {
    auto s = state();
    auto new_name = name_+"["+std::to_string(n)+"]";
    // push the table onto the stack
    push_last();
    // make sure we are accessing a table
    check_table(name_);
    // now get the table value, the key gets pushed from the stack
    lua_rawgeti(s, -1, n);
    // make sure returned result is not nil
    check_result(new_name);
    // get a reference to the result, and pop the table from the stack
    auto ref = make_lua_ref(state_);
    lua_pop(s, -1);
    // return the global object with a pointer to a location in the stack
    return { state_, new_name, std::move(ref) };
  }

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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
    return { state_, name, make_lua_ref(state_) };
  }

  auto operator()( const std::string & script )
  {
    return run_string( script );
  }
};

#endif
  
} // namespace utils
} // namespace ale
