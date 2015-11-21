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
 * \file check_types.h
 * 
 * \brief Some helper functions to check variable types.
 *
 ******************************************************************************/
#pragma once

namespace ale {
namespace utils {


////////////////////////////////////////////////////////////////////////////////
//! \brief Test to see if all variadic template arguments are of type T
////////////////////////////////////////////////////////////////////////////////

template<typename... Conds>
struct and_ : std::true_type {};

// combine conditions
template<typename Cond, typename... Conds>
struct and_<Cond, Conds...> : 
  std::conditional< Cond::value, 
                    and_<Conds...>,
                    std::false_type >::type 
{};

// This is the exposed function!!
template<typename Target, typename... Ts>
using are_type_t = and_<std::is_same<Ts,Target>...>;


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
