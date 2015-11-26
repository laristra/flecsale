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
 * \file tuple_visit.h
 * 
 * \brief Some helper functions for looping through tuples.
 *
 ******************************************************************************/
#pragma once

namespace ale {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! \brief zip or tie two tuples together
////////////////////////////////////////////////////////////////////////////////

// support struct to iterate over the tuple(s)
template<size_t size>
struct visit_tuple_ws
{
  template<typename Callable, typename Head, typename... Tail>
  static void visit(Callable&& f, Head&& aTuple, Tail&&... aTail)
  {
    visit_tuple_ws<size-1>::visit( std::forward<Callable>(f), 
                                   std::forward<Head>(aTuple), 
                                   std::forward<Tail>(aTail)... );
    f( std::get<size>( std::forward<Head>(aTuple) ), 
       std::get<size>( std::forward<Tail>(aTail) )...);
  }
};

// stop recursion here
template<>
struct visit_tuple_ws<0u>
{
  template<typename Callable, typename Head, typename... Tail>
  static void visit(Callable&& f, Head&& aTuple, Tail&&... aTail)
  {
    f( std::get<0>( std::forward<Head>(aTuple) ), 
       std::get<0>( std::forward<Tail>(aTail) )...);
  }
};

// visit_tuple
template<typename Callable, typename Head, typename... Tail>
Callable tuple_visit(Callable&& f, Head&& aTuple, Tail&&... aTail)
{
  const size_t size = std::tuple_size<typename std::remove_reference<Head>::type>::value;
  visit_tuple_ws<size-1>::visit( std::forward<Callable>(f), 
                                 std::forward<Head>(aTuple), 
                                 std::forward<Tail>(aTail)... );
  return f;
}


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
