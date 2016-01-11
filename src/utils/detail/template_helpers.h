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

namespace ale {
namespace utils {

namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Template helper to statically multiply arguments together
////////////////////////////////////////////////////////////////////////////////


//! \brief return 1 for the final multiplcation
constexpr std::size_t multiply() 
{ return 1; }

//! \brief main implementation for multiplication
template<typename Arg, typename... Args>
constexpr auto multiply(Arg first, Args... rest) 
{ return first * multiply(rest...); }


}  // namespace


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
