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
//! \remark this is my version
////////////////////////////////////////////////////////////////////////////////

#if 1


// create index list
template<class T>
constexpr std::make_index_sequence<std::tuple_size<T>::value>
get_indexes( T const& )
{ return {}; }

// actuall call to functions
template<size_t... Is, class Tuple, class F>
void tuple_for_each( std::index_sequence<Is...>, Tuple&& tup, F&& f ) {
  using std::get;
  int unused[] = { 0, ( (void)f(get<Is>(std::forward<Tuple>(tup))), 0 )... };
}

// This is the exposed function!!
template<class Tuple, class F>
void tuple_for_each( Tuple&& tup, F&& f ) {
  auto indexes = get_indexes(tup);
  tuple_for_each(indexes, std::forward<Tuple>(tup), std::forward<F>(f) );
}


#else

////////////////////////////////////////////////////////////////////////////////
//! \brief Exectute something for each element of a tuple
//! \remark this is ben's version
////////////////////////////////////////////////////////////////////////////////

// convenience define
template<typename TupleType>
using __flexi_integral_constant_t = std::integral_constant<size_t,
   std::tuple_size<typename std::remove_reference<TupleType>::type>::value>;

// Empty sequence call to terminate
template<typename TupleType, typename FunctionType>
void tuple_for_each(TupleType &&, FunctionType,
   __flexi_integral_constant_t<TupleType>) {}

// Expand each index
template<
   size_t I,
   typename TupleType,
   typename FunctionType,
   typename = typename std::enable_if<I!=std::tuple_size<
      typename std::remove_reference<TupleType>::type>::value>::type
>
void tuple_for_each(TupleType && t, FunctionType f,
   std::integral_constant<size_t, I>) {
      f(std::get<I>(t));
      tuple_for_each(std::forward<TupleType>(t), f,
         std::integral_constant<size_t, I + 1>());
}

// Main interface
template<typename TupleType, typename FunctionType>
void tuple_for_each(TupleType && t, FunctionType f)
{
    tuple_for_each(std::forward<TupleType>(t), f,
      std::integral_constant<size_t, 0>());
}

#endif


} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/

