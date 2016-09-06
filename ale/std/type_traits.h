/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Extensions to C++11's <type_traits> library that are part of C++14.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <type_traits>

namespace std {

////////////////////////////////////////////////////////////////////////////////
// Stuff not available till C++14
////////////////////////////////////////////////////////////////////////////////

#if __cplusplus == 201103L

//==============================================================================
//! \brief An enable_if helper type.
//! \tparam B  A boolean operation, if true, convert to type T.
//! \tparam T  The type to convert to upon a true boolean operation.
//==============================================================================
template< bool B, class T = void >
using enable_if_t = typename enable_if<B,T>::type;

//==============================================================================
//! \brief A remove_reference helper type.
//! \tparam T  The type to remove references from.
//==============================================================================
template< class T >
using remove_reference_t = typename remove_reference<T>::type;

//==============================================================================
//! \brief A remove_const helper type.
//! \tparam T  The type to remove const from.
//==============================================================================
template< class T >
using remove_const_t = typename remove_const<T>::type;


//==============================================================================
//! \brief a remove_pointer helper type.
//! \tparam T  The type to remove the pointer from.
//==============================================================================
template< class T >
using remove_pointer_t = typename remove_pointer<T>::type;

//==============================================================================
//! \brief A decay helper type.
//! \tparam T  The type to decay.
//==============================================================================
template< class T >
using decay_t = typename decay<T>::type;

//==============================================================================
//! \brief A is_pointer helper type.  True if T is a pointer type.
//! \tparam T  The type to test if it is a pointer.
//==============================================================================
template< class T >
constexpr bool is_pointer_v = is_pointer<T>::value;

////////////////////////////////////////////////////////////////////////////////
// Stuff not available till C++17
////////////////////////////////////////////////////////////////////////////////

#elif __cplusplus == 201402L

//==============================================================================
//! \brief A is_same helper type.  True if T and U are the same types.
//! \tparam T,U  The types to test if they are the same.
//==============================================================================
template< class T, class U >
constexpr bool is_same_v = is_same<T, U>::value;


//==============================================================================
//! \brief A is_arithmetic helper type.  True if T is an arithmetic type.
//! \tparam T  The type to test.
//==============================================================================
template< class T >
constexpr bool is_arithmetic_v = is_arithmetic<T>::value;

//==============================================================================
//! \brief A is_function helper type.  True if T is a std::function.
//! \tparam T  The type to test.
//==============================================================================
template< class T >
constexpr bool is_function_v = is_function<T>::value;


//==============================================================================
//! \brief A is_integral helper type.  True if T is an integral type.
//! \tparam T  The type to test.
//==============================================================================
template< class T >
constexpr bool is_integral_v = is_integral<T>::value;

#endif

} // namespace
