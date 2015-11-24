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
#include "ale/utils/static_for_each.h"
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
fill( const tuple_t<Types...> t, const auto& val ) 
{ 
  //std::cout << "tuple_t (single value constructor)\n";
  utils::static_for_each( t, [&](auto & tup) { tup = val; } );
}


//! \brief Equivalence operator
//! \param[in] lhs The quantity on the rhs.
//! \param[in] rhs The quantity on the rhs.
//! \return true if equality.
template <typename... Types>
bool operator==(const tuple_t& lhs, const tuple_t& rhs)
{
  bool result = true;
  if ( &lhs != &rhs ) {
    utils::static_for_each( utils::tuple_tie( lhs.data_, rhs.data_ ),
                     [&](auto && tup) { 
                       if ( std::get<0>(tup) != std::get<1>(tup) )
                         result = false;
                     } );    
    return result;
  }
}



//! \brief Output operator for tuple_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in,out] os  The ostream to dump output to.
//! \param[in]     rhs The tuple_t on the right hand side of the operator.
//! \return A reference to the current ostream.
template <typename... Types>
auto & operator<<(std::ostream& os, const tuple_t& a)
{
  os << "{";
  utils::static_for_each( a.data_, 
                   [](auto & tup) { std::cout << " [ " << tup << " ]"; } );
  os << " }";
  return os;
}
  
 
 
//! \brief Addition operator involving two tuple_ts.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The tuple_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator+( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp += rhs;
  return tmp;
}


//! \brief Addition operator involving one tuple_t and a scalar.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator+( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp += rhs;
  return tmp;
}

template <typename... Types>
auto operator+( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp += rhs;
  return tmp;
}

//! \brief Subtraction operator involving two tuple_ts.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The tuple_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator-( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

//! \brief Subtraction operator involving one tuple_t and a scalar.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator-( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

template <typename... Types>
auto operator-( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

//! \brief Multiplication operator involving two tuple_ts.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The tuple_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator*( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp *= rhs;
  return tmp;
}


//! \brief Multiplication operator involving one tuple_t and a scalar.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator*( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

template <typename... Types>
auto operator*( const auto & lhs,
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

//! \brief Division operator involving two tuple_ts.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The tuple_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator/( const tuple_t<Types...>& lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp /= rhs;
  return tmp;
}



//! \brief Division operator involving one tuple_t and a scalar.
//! \param[in] lhs The tuple_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
auto operator/( const tuple_t<Types...>& lhs, 
                const auto & rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp /= rhs;
  return tmp;
}

template <typename... Types>
auto operator/( const auto & lhs, 
                const tuple_t<Types...>& rhs )
{
  tuple_t<Types...> tmp(lhs);
  tmp /= rhs;
  return tmp;
}








#if 0
  
//! \brief Addition binary operator involving another array.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename... Types>
void add( const T&& lhs, const tuple_t<Types...> &rhs) {
  utils::static_for_each( utils::tuple_tie( lhs, rhs ),
                          [](auto && tup) { 
                            std::get<0>(tup) += std::get<1>(tup);
                          } );
}
  

//! \brief Multiplication binary operator involving a constant.
//! \param[in] val The constant on the right hand side of the operator.
//! \return A reference to the current object.
tuple_t & operator*=(const auto &val) {
  utils::static_for_each( data_,
                   [&](auto && tup) { tup *= val; } );
  return *this;
}

//! \brief Division binary operator involving another array.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
auto & operator/=(const tuple_t &rhs) {
  if ( this != &rhs ) {
    utils::static_for_each( utils::tuple_tie( data_, rhs.data_ ),
                     [](auto && tup) { 
                       std::get<0>(tup) /= std::get<1>(tup);
                     } );
  }
  return *this;
}

//! \brief Division operator involving a constant.
//! \param[in] val The constant on the right hand side of the operator.
//! \return A reference to the current object.
tuple_t & operator/=(const auto &val) {
  utils::static_for_each( data_,
                   [&](auto && tup) { tup /= val; } );
  return *this;
}




//! \brief index operator
template<size_t I, class... UTypes>
auto & get ( tuple_t<UTypes...> & tup ) 
{ return std::get<I>(tup.data_); }
  

//! \brief index operator (const version)
template<size_t I, class... UTypes>
const auto & get ( const tuple_t<UTypes...> & tup ) 
{ return std::get<I>(tup.data_); }


//! \brief Extract the tuple element type.
//! \remark undefined
template <size_t I, class T>
struct tuple_element;

//! \brief Extract the tuple element type.
template <size_t I, class... Types>
struct tuple_element<I, tuple_t<Types...>> 
{
private:
  using data_t = typename tuple_t<Types...>::data_t;
public:
  using type = typename std::tuple_element<I, data_t>::type;
};

//! \brief Extract the tuple element size.
//! \remark undefined
template <class T>
struct tuple_size;

//! \brief Extract the tuple size.
template <class... Types>
struct tuple_size<tuple_t<Types...>> 
{
private:
  using data_t = typename tuple_t<Types...>::data_t;
public:
  using value = typename std::tuple_size<data_t>::value;
};



//! \brief Copy tuple's forward_as_tuple functionality.
template<class... Types>
tuple_t<Types&&...> forward_as_tuple (Types&&... args) noexcept
{
  return tuple_t<Types&&...>(std::forward<Types>(args)...);
}

#endif

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
