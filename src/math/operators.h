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

//! \brief Fill a container with one value.
//! \param[in] val The value to set the array to
template< class T = void >
struct fill;


template< class T >
struct fill {
  constexpr void operator()( T & lhs, const T & rhs ) 
  { lhs = rhs; }
};


//! \brief Addition operator.
template< class T = void >
struct plus;

template<>
struct plus<void> {
  template< class T, class U >
  constexpr auto operator()( T && lhs, U && rhs ) 
  { return std::forward<T>(lhs) + std::forward<U>(rhs); }
};

//! \brief Addition assigment operator.
template< class T = void >
struct plus_equal;

template<>
struct plus_equal<void> {
  template< class T, class U >
  constexpr void operator()( T && lhs, T && rhs ) 
  { std::forward<T>(lhs) += std::forward<U>(rhs); }
};

//! \brief Minus operator.
template< class T = void >
struct minus;

template< class T >
struct minus {
  constexpr auto operator()( const T & lhs, const T & rhs ) 
  { return lhs - rhs; }
};
 
//! \brief Minus assigment operator.
template< class T = void >
struct minus_equal;

template< class T >
struct minus_equal {
  constexpr void operator()( T & lhs, const T & rhs ) 
  { lhs -= rhs; }
};

//! \brief Multiplication operator.
template< class T = void >
struct multiplies;

template< class T >
struct multiplies {
  constexpr auto operator()( const T & lhs, const T & rhs ) 
  { return lhs * rhs; }
};

//! \brief Multiplication assigment operator.
template< class T = void >
struct multiplies_equal;

template< class T >
struct multiplies_equal {
  constexpr void operator()( T & lhs, const T & rhs ) 
  { lhs *= rhs; }
};


//! \brief Division operator.
template< class T = void >
struct divides;

template< class T >
struct divides {
  constexpr auto operator()( const T & lhs, const T & rhs ) 
  { return lhs / rhs; }
};

//! \brief Division assigment operator.
template< class T = void >
struct divides_equal;

template< class T >
struct divides_equal {
  constexpr void operator()( T & lhs, const T & rhs ) 
  { lhs /= rhs; }
};


} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
