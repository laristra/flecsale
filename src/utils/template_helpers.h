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
