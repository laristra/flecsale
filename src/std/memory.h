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
 * \file memory.h
 * 
 * \brief Extensions to C++11 that are part of C++14.
 *
 ******************************************************************************/
#pragma once

//! system includes
#include <memory>

#if __cplusplus == 201103L

//! system includes needed for this implementation
#include <cstddef>
#include <type_traits>
#include <utility>

////////////////////////////////////////////////////////////////////////////////
//! \brief an implementation of make_unique
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
////////////////////////////////////////////////////////////////////////////////


namespace std {

  template<class T> struct _Unique_if {
    typedef unique_ptr<T> _Single_object;
  };

  template<class T> struct _Unique_if<T[]> {
    typedef unique_ptr<T[]> _Unknown_bound;
  };

  template<class T, size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
  };

  //! \brief Constructs a non-array type T
  //! \param [in]  args  the arguments to the constructor
  //! \return a unique_ptr
  template<class T, class... Args>
  typename _Unique_if<T>::_Single_object
  make_unique(Args&&... args) {
    return unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  //! \brief Constructs an array of unknown bound T.
  //! \param [in]  n  the size of the array
  //! \return a unique_ptr
  template<class T>
  typename _Unique_if<T>::_Unknown_bound
  make_unique(size_t n) {
    typedef typename remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
  }

  //! \brief Construction of arrays of known bound is disallowed.
  template<class T, class... Args>
  typename _Unique_if<T>::_Known_bound
  make_unique(Args&&...) = delete;

} // namespace std


#endif

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
