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
 * \file functional.h
 * 
 * \brief Extends some functionality of the Functional library.
 *
 ******************************************************************************/
#pragma once

//! system includes
#include <functional>

namespace std {

////////////////////////////////////////////////////////////////////////////////
// Allow functional operators to work for two different argument types
////////////////////////////////////////////////////////////////////////////////

//! \brief Addition operator.
template< class T, class U >
constexpr auto plus( T && lhs, U && rhs )
{ return std::forward<T>(lhs) + std::forward<U>(rhs); }

//! \brief Minus operator.
template< class T, class U >
constexpr auto minus( T && lhs, U && rhs )
{ return std::forward<T>(lhs) - std::forward<U>(rhs); }
 
//! \brief Multiplication operator.
template< class T, class U >
constexpr auto multiplies( T && lhs, U && rhs )
{ return std::forward<T>(lhs) * std::forward<U>(rhs); }

//! \brief Division operator.
template< class T, class U >
constexpr auto divides( T && lhs, U && rhs )
{ return std::forward<T>(lhs) / std::forward<U>(rhs); }


} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
