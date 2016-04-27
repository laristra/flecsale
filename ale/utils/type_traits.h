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
 * \brief Statically check if all arguments are of the same type.
 *
 ******************************************************************************/
#pragma once

// system inculdes
#include <type_traits>

// user includes
#include "detail/type_traits.h"


namespace ale {
namespace utils {


////////////////////////////////////////////////////////////////////////////////
//! \brief Test to see if all variadic template arguments are of type Target
////////////////////////////////////////////////////////////////////////////////
template<typename Target, typename... Ts>
using are_type_t = detail::and_< 
  std::is_same< std::decay_t<Ts>, std::decay_t<Target> >... >;



////////////////////////////////////////////////////////////////////////////////
//! \brief check if a particular type is an iterator
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct is_iterator {
  static char test(...);

  template <typename U,
            typename=typename std::iterator_traits<U>::difference_type,
            typename=typename std::iterator_traits<U>::pointer,
            typename=typename std::iterator_traits<U>::reference,
            typename=typename std::iterator_traits<U>::value_type,
            typename=typename std::iterator_traits<U>::iterator_category
  > static long test(U&&);

  constexpr static bool value = std::is_same<decltype(test(std::declval<T>())),long>::value;
};

////////////////////////////////////////////////////////////////////////////////
//! \brief helper function for checking if iterator
////////////////////////////////////////////////////////////////////////////////
template< typename T >
constexpr bool is_iterator_v = is_iterator<T>::value;

} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
