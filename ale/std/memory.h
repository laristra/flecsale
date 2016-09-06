/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Extensions to C++11's <memory> library  that are part of C++14.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <memory>

////////////////////////////////////////////////////////////////////////////////
// Stuff not available till C++14
////////////////////////////////////////////////////////////////////////////////

#if __cplusplus == 201103L

// system includes needed for this implementation
#include <cstddef>
#include <type_traits>
#include <utility>


namespace std {

namespace detail {

//! \brief Instantiated if T is a single object.
template<class T> struct _Unique_if {
  typedef unique_ptr<T> _Single_object;
};

//! \brief Instantiated if T is a array with unknown bounds.
template<class T> struct _Unique_if<T[]> {
  typedef unique_ptr<T[]> _Unknown_bound;
};

//! \brief Instantiated if T is an array with known bounds.
template<class T, size_t N> struct _Unique_if<T[N]> {
  typedef void _Known_bound;
};

} // namespace detail

////////////////////////////////////////////////////////////////////////////////
//! \defgroup make_unique Various implementations of std::make_unique.
//! 
//! \brief An implementation of make_unique.
//!
//! Constructs an object of type T and wraps it in a std::unique_ptr. 
//!
//! 1. Constructs a non-array type T. The arguments args are passed to the 
//!    constructor of T. The function does not participate in the overload 
//!    resolution if T is an array type. The function is equivalent to:
//!      unique_ptr<T>(new T(std::forward<Args>(args)...))
//! 
//! 2. Constructs an array of unknown bound T. The function does not participate 
//!    in the overload resolution unless T is an array of unknown bound. The 
//!    function is equivalent to: 
//!      unique_ptr<T>(new typename std::remove_extent<T>::type[size]())
//!
//! 3. Construction of arrays of known bound is disallowed.
//!
//! \remark Not needed in C++14
//! \return Return a unique_ptr to the object.
//! \tparam T The non-array type to construct.
////////////////////////////////////////////////////////////////////////////////
///@{ 

//! \brief This version constructs a non-array type T.
//! \param [in]  args  The arguments to the constructor.
//! \tparam Args The variadic argument types passed to the objects constructor.
template<class T, class... Args>
typename detail::_Unique_if<T>::_Single_object
make_unique(Args&&... args) {
  return unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//! \brief This version constructs an array of unknown bound T.
//! \param [in]  n  The size of the array.
template<class T>
typename detail::_Unique_if<T>::_Unknown_bound
make_unique(size_t n) {
  typedef typename remove_extent<T>::type U;
  return unique_ptr<T>(new U[n]());
}

//! \brief  Construction of arrays of known bound is disallowed.
//! \tparam Args The variadic argument types passed to the objects constructor.
template<class T, class... Args>
typename detail::_Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;

///@}
} // namespace std


#endif
