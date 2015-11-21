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
 * \file tuple_zip.h
 * 
 * \brief Some helper functions for zipping tupples together.
 *
 ******************************************************************************/
#pragma once

namespace ale {
namespace utils {


////////////////////////////////////////////////////////////////////////////////
//! \brief zip or tie two tuples together
////////////////////////////////////////////////////////////////////////////////

// helpers for generating sequences
template<size_t ...S>
struct seq { };

template<size_t N, size_t ...S>
struct gens : gens<N-1, N-1, S...> 
{ };

template<size_t ...S>
struct gens<0, S...> {
  typedef seq<S...> type;
};

// zip each tuple together
template <class Tup1, class Tup2, size_t ...S>
auto tuple_zip_helper(Tup1&& t1, Tup2&& t2, seq<S...> s)  {
  return std::make_tuple( std::make_pair( std::get<S>(t1), std::get<S>(t2) )...);
}



// main interface zip function using perfect forwarding
template <class Tup1, class Tup2>
auto tuple_zip(Tup1&& t1, Tup2&& t2)  {
  constexpr auto num_args1 = std::tuple_size<std::decay_t<Tup1>>::value;
  constexpr auto num_args2 = std::tuple_size<std::decay_t<Tup2>>::value;
  static_assert(num_args1 == num_args2, "The tuple sizes must be the same");
  return tuple_zip_helper( std::forward<Tup1>(t1), std::forward<Tup2>(t2), 
                           typename gens<num_args2>::type() );
}

// tie each tuple together
template <class Tup1, class Tup2, size_t ...S>
auto tuple_tie_helper(Tup1&& t1, Tup2&& t2, seq<S...> s)  {
  return std::make_tuple( std::forward_as_tuple( std::get<S>(t1), std::get<S>(t2) )...);
}



// main interface tie function using perfect forwarding
template <class Tup1, class Tup2>
auto tuple_tie(Tup1&& t1, Tup2&& t2)  {
  constexpr auto num_args1 = std::tuple_size<std::decay_t<Tup1>>::value;
  constexpr auto num_args2 = std::tuple_size<std::decay_t<Tup2>>::value;
  static_assert(num_args1 == num_args2, "The tuple sizes must be the same");
  return tuple_tie_helper( std::forward<Tup1>(t1), std::forward<Tup2>(t2), 
                           typename gens<num_args2>::type() );
}



} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
