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

// user includes
#include "detail/template_helpers.h"

namespace ale {
namespace utils {



////////////////////////////////////////////////////////////////////////////////
//! \brief statically multiply arguments together
////////////////////////////////////////////////////////////////////////////////
template<typename... Args>
constexpr auto multiply(Args... args) 
{ return detail::multiply(args...); }


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
