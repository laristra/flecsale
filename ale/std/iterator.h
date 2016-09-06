/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Extends some functionality of the <iterator> library.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <iterator>

namespace std {

////////////////////////////////////////////////////////////////////////////////
// Allow begin and end to work for scalar values
////////////////////////////////////////////////////////////////////////////////

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

} // namespace
