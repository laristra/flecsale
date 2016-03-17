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
 * \file const_str.h
 * 
 * \brief A constexpr string.
 *
 * \see http://en.cppreference.com/w/cpp/language/constexpr
 * \see https://github.com/boostcon/cppnow_presentations_2012/blob/master/wed/schurr_cpp11_tools_for_class_authors.pdf?raw=true
 ******************************************************************************/
#pragma once

namespace ale {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! \brief constexpr string
////////////////////////////////////////////////////////////////////////////////
class const_string {

public:

  using size_type = std::size_t;
  using hash_type = size_type;

  //! \brief constructor
  template< size_type N >
  constexpr const_string( const char (&a)[N] ) : p_(a), sz_(N-1) {}

  //! \brief operator []
  constexpr char operator[] ( size_type n )
  { return n < sz_ ? p_[n] : throw std::out_of_range(""); }

  //! \brief size()
  constexpr size_type size() { return sz_; }

  //! \brief c_str accessor
  constexpr const char* c_str() const { return p_; }

private:

  const char* const p_;
  const size_type sz_;

};

} // namspeace
} // namspeace
