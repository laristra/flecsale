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
 * \file operators.h
 * 
 * \brief Provides a default operators for fundamental types.
 *
 ******************************************************************************/
#pragma once


namespace ale {
namespace math {


//! \brief User-defined overloaded begin/end for scalar values
//! \param[in] a the scalar
//! \return the pointer to the beginning
template <typename T>
typename std::enable_if< std::is_scalar<T>::value, T* >::type 
begin( T & a ) { return &a; };

//! \brief User-defined overloaded begin/end for scalar values
//! \param[in] a the scalar
//! \return the pointer to the end
template <typename T>
typename std::enable_if< std::is_scalar<T>::value, T* >::type 
end( T & a ) { return &a+1; };


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
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
