/*~-------------------------------------------------------------------------~~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
/*!
 *
 * \file template_helpers.h
 * 
 * \brief Some helper functions for template foo magic.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <functional>

namespace ale {
namespace utils {



////////////////////////////////////////////////////////////////////////////////
//! \brief statically multiply arguments together
////////////////////////////////////////////////////////////////////////////////

namespace detail {

//! \brief return 1 for the final multiplcation
constexpr std::size_t multiply() 
{ return 1; }

//! \brief main implementation for multiplication
template<typename Arg, typename... Args>
constexpr auto multiply(Arg first, Args... rest) 
{ return first * multiply(rest...); }


}  // namespace

//! \brief the main interface
template<typename... Args>
constexpr auto multiply(Args... args) 
{ return detail::multiply(args...); }

////////////////////////////////////////////////////////////////////////////////
//! \brief statically fill an array with a constant value
////////////////////////////////////////////////////////////////////////////////

template <std::size_t N>
struct fill {
  template <typename T, typename ...Tn>
  static constexpr auto apply(T v, Tn ...vs)
  {
    return fill<N - 1>::apply(v, v, vs...);
  }
};

template <>
struct fill<1> {
  template <typename T, typename ...Tn>
  static constexpr auto apply(T v, Tn ...vs)
  {
    return std::array<T, sizeof...(vs) + 1>{v, vs...};
  }

};

////////////////////////////////////////////////////////////////////////////////
//! \brief statically make an array with a constant value
////////////////////////////////////////////////////////////////////////////////

namespace detail {

template <typename T, std::size_t...Is>
constexpr std::array<T, sizeof...(Is)> make_array(T val, std::index_sequence<Is...>)
{
  return {(static_cast<void>(Is), val)...};
}

} // namespace 


template <typename T, std::size_t N>
constexpr std::array<T, N> make_array(T val)
{
  return detail::make_array(val, std::make_index_sequence<N>());
}

////////////////////////////////////////////////////////////////////////////////
//! \brief a tie using constant references
////////////////////////////////////////////////////////////////////////////////
template < typename T, typename... Ts >
std::tuple<T&, const Ts&...> ctie( T& first, const Ts&... rest )
{
  return std::make_tuple( std::ref(first), std::cref(rest)... );
}

////////////////////////////////////////////////////////////////////////////////
// A type traits type struct used to unpack a tuple type, and repack it using 
// references
////////////////////////////////////////////////////////////////////////////////

//! \brief Unpack a tuple and create a tuple of references to each element.
//! \tparam T  The tuple type
//! \remark This is the empty struct
template < typename T >
struct reference_wrapper {};

//! \brief Unpack a tuple and create a tuple of references to each element.
//! \tparam Tuple  The tuple type
//! \remark This is the tuple implementation
template < template<typename...> class Tuple, typename... Args >
struct reference_wrapper < Tuple<Args...> >
{
  using type = Tuple<Args&...>;
};


//! \brief This is a helper function to unpack a tuple and create references
//! \tparam T  The tuple type
template < typename T >
using reference_wrapper_t = typename reference_wrapper<T>::type;


////////////////////////////////////////////////////////////////////////////////
//! \brief return an lvalue reference
////////////////////////////////////////////////////////////////////////////////
template<typename T>
T &as_lvalue(T &&val) {
  return val;
}

} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
