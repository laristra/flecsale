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
 * \file tuple.h
 * 
 * \brief Provides a tuple type which functions as a vector.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <tuple>

//! user includes
#include "ale/utils/tuple_for_each.h"
#include "ale/utils/tuple_visit.h"
#include "ale/utils/tuple_zip.h"

namespace ale {
namespace math {

////////////////////////////////////////////////////////////////////////////////
//!  \brief The dimensioned_array type provides a general base for defining
//!  contiguous array types that have a specific dimension.
//!
//!  \tparam T The type of the array, e.g., P.O.D. type.
//!  \tparam D The dimension of the array, i.e., the number of elements
//!    to be stored in the array.
////////////////////////////////////////////////////////////////////////////////
template <typename... Types>
using tuple_t = std::tuple<Types...>;

//! \brief Constructor with one value.
//! \param[in] val The value to set the array to
template < typename... Types >
void fill( tuple_t<Types...> & t, const auto& val ) 
{ 
  //std::cout << "tuple_t (single value constructor)\n";
  utils::tuple_for_each( t, [&](auto & tup) { tup = val; } );
}


//! \brief Add to a tuple_t.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename... Types>
void add_to( tuple_t<Types...>& lhs, 
             const tuple_t<Types...>& rhs )
{
  utils::tuple_for_each( utils::tuple_tie( lhs, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<0>(tup) + std::get<1>(tup);
                         } );                        
}

template <typename... Types>
void add_to( tuple_t<Types...>& lhs, 
             const auto & rhs )
{
  utils::tuple_for_each( lhs,
                         [&](auto & tup) { 
                           tup = tup + rhs;
                         } );                        
}


//! \brief Addition operator involving tuple_ts.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename... Types>
auto operator+( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_visit( 
                     [](auto & a, const auto & b, const auto & c) { 
                       a = b + c;
                     }, 
                     tmp, lhs, rhs );
  return tmp;
}

template <typename... Types>
auto operator+( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, lhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<1>(tup) + rhs;
                         } );                        
  return tmp;
}

template <typename... Types>
auto operator+( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = lhs + std::get<1>(tup);
                         } );                        
  return tmp;
}



//! \brief Subtract a tuple_t.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename... Types>
void subtract_from( tuple_t<Types...>& lhs, 
                    const tuple_t<Types...>& rhs )
{
  utils::tuple_for_each( utils::tuple_tie( lhs, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<0>(tup) - std::get<1>(tup);
                         } );                        
}

template <typename... Types>
void subtract_from( tuple_t<Types...>& lhs, 
                    const auto & rhs )
{
  utils::tuple_for_each( lhs,
                         [&](auto & tup) { 
                           tup = tup - rhs;
                         } );                        
}


//! \brief Subtraction operator involving tuple_ts.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename... Types>
auto operator-( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_visit( 
                     [](auto & a, const auto & b, const auto & c) { 
                       a = b - c;
                     }, 
                     tmp, lhs, rhs );
  return tmp;
}

template <typename... Types>
auto operator-( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, lhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<1>(tup) - rhs;
                         } );                        
  return tmp;
}

template <typename... Types>
auto operator-( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = lhs - std::get<1>(tup);
                         } );                        
  return tmp;
}


//! \brief Multiply a tuple_t.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename... Types>
void multiply_by( tuple_t<Types...>& lhs, 
                  const tuple_t<Types...>& rhs )
{
  utils::tuple_for_each( utils::tuple_tie( lhs, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<0>(tup) * std::get<1>(tup);
                         } );                        
}

template <typename... Types>
void multiply_by( tuple_t<Types...>& lhs, 
                  const auto & rhs )
{
  utils::tuple_for_each( lhs,
                         [&](auto & tup) { 
                           tup = tup * rhs;
                         } );                        
}


//! \brief Multiplication operator involving tuple_ts.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename... Types>
auto operator*( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_visit( 
                     [](auto & a, const auto & b, const auto & c) { 
                       a = b * c;
                     }, 
                     tmp, lhs, rhs );
  return tmp;
}

template <typename... Types>
auto operator*( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, lhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<1>(tup) * rhs;
                         } );                        
  return tmp;
}

template <typename... Types>
auto operator*( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = lhs * std::get<1>(tup);
                         } );                        
  return tmp;
}





//! \brief Divide a tuple_t.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
template <typename... Types>
void divide_by( tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  utils::tuple_for_each( utils::tuple_tie( lhs, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<0>(tup) / std::get<1>(tup);
                         } );                        
}

template <typename... Types>
void divide_by( tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  utils::tuple_for_each( lhs,
                         [&](auto & tup) { 
                           tup = tup / rhs;
                         } );                        
}


//! \brief Division operator involving tuple_ts.
//! \param[in] lhs The value on the left hand side of the operator.
//! \param[in] rhs The value on the right hand side of the operator.
//! \return The result of the operation.
template <typename... Types>
auto operator/( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_visit( 
                     [](auto & a, const auto & b, const auto & c) { 
                       a = b / c;
                     }, 
                     tmp, lhs, rhs );
  return tmp;
}

template <typename... Types>
auto operator/( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, lhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = std::get<1>(tup) / rhs;
                         } );                        
  return tmp;
}

template <typename... Types>
auto operator/( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp;
  utils::tuple_for_each( utils::tuple_tie( tmp, rhs ),
                         [&](auto tup) { 
                           std::get<0>(tup) = lhs / std::get<1>(tup);
                         } );                        
  return tmp;
}




//! \brief Output operator for tuple_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in,out] os  The ostream to dump output to.
//! \param[in]     rhs The tuple_t on the right hand side of the operator.
//! \return A reference to the current ostream.
template <typename... Types>
auto & operator<<(std::ostream& os, const tuple_t<Types...>& a)
{
  os << "{";
  utils::tuple_for_each( a, 
                          [](auto & tup) { std::cout << " [ " << tup << " ]"; } );
  os << " }";
  return os;
}
  
 
 

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
