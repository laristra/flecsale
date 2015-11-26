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
 * \file tuple_for_each.h
 * 
 * \brief A static for-each function for looping over tuples statically.
 *
 ******************************************************************************/
#pragma once

namespace ale {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! \brief Exectute something for each element of a tuple
//! \remark this is ben's version
////////////////////////////////////////////////////////////////////////////////


// actuall call to functions
template<size_t... Is, class F>
void static_for( std::index_sequence<Is...>, F&& f ) {
  int unused[] = { 0, ( (void)f(Is), 0 )... };
}

// This is the exposed function!!
template<size_t N, class F>
void static_for(  F&& f ) {
  auto indexes = std::make_index_sequence<N>();
  static_for(indexes, std::forward<F>(f) );
}



} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/

