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
template <typename... Types> class tuple_t {

public:

  //===========================================================================
  // Typedefs
  //===========================================================================

  //! \brief the size_t type
  using size_t = std::size_t;

  //! the size of the tupel
  static constexpr size_t size = sizeof...(Types);

  //! alias the tuple data type
  using tuple_data_t = std::tuple<Types...>;

  //===========================================================================
  // Constructors / Destructors
  //===========================================================================

  //! \brief force the default constructor
  tuple_t() = default;

  //! \brief force the default copy constructor
  tuple_t(const tuple_t &) = default;


  //! \brief Constructor with initializer list
  //! \param[in] list the initializer list of values
  template < typename... UTypes,
             typename = typename std::enable_if< 
               sizeof...(UTypes) == size &&
               sizeof...(UTypes) >= 2 >::type
             >
  tuple_t(UTypes&&... args) : data_(std::forward<UTypes>(args)...)
  { 
    //std::cout << "tuple_t (variadic constructor)\n";
  }


  //! \brief Constructor with one value.
  //! \param[in] val The value to set the array to
  tuple_t(const auto & val) 
  { 
    //std::cout << "tuple_t (single value constructor)\n";
    utils::static_for_each( data_, [&](auto & tup) { tup = val; } );
  }


  //===========================================================================
  // Operators
  //===========================================================================

  //! \brief index operator
  template<size_t I>
  auto & get()
  { return std::get<I>(data_); }

  //! \brief index operator (forwarding version)
  template<size_t I, class... ArgTypes>
  friend auto && get( tuple_t<ArgTypes...> && tup );

  //! \brief index operator (const version)
  template<size_t I>
  const auto & get() const 
  {  return std::get<I>(data_); }


  //! \brief Assignment operator to another array.
  //! \param[in] rhs The array on the right hand side of the '='.
  //! \return A reference to the current object.
  auto & operator=(const tuple_t &rhs) {
    if ( this != &rhs ) data_ = rhs.data_;
    return *this;
  }

  //! \brief Assignment operator to a constant.
  //! \param[in] val The constant on the right hand side of the '='.
  //! \return A reference to the current object.
  tuple_t & operator=(const auto & val) {
    utils::static_for_each( data_, [&](auto & tup) { tup = val; } );
    return *this;
  }
  
  //! \brief Addition binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator+=(const tuple_t &rhs) {
    if ( this != &rhs ) {
      utils::static_for_each( utils::tuple_tie( data_, rhs.data_ ),
                       [](auto && tup) { 
                         std::get<0>(tup) += std::get<1>(tup);
                       } );
    }
    return *this;
  }
  
  //! \brief Addiition binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  tuple_t & operator+=(const auto &val) {
    utils::static_for_each( data_,
                     [&](auto && tup) { tup += val; } );
    return *this;
  }

  //! \brief Subtraction binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator-=(const tuple_t &rhs) {
    if ( this != &rhs ) {
      utils::static_for_each( utils::tuple_tie( data_, rhs.data_ ),
                       [](auto && tup) { 
                         std::get<0>(tup) -= std::get<1>(tup);
                       } );
    }
    return *this;
  }

  //! \brief Subtraction binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  tuple_t & operator-=(const auto &val) {
    utils::static_for_each( data_,
                     [&](auto && tup) { tup -= val; } );
    return *this;
  }


  //! \brief Multiplication binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator*=(const tuple_t &rhs) {
    if ( this != &rhs ) {
      utils::static_for_each( utils::tuple_tie( data_, rhs.data_ ),
                       [](auto && tup) { 
                         std::get<0>(tup) *= std::get<1>(tup);
                       } );
    }
    return *this;
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

  //! \brief Equivalence operator
  //! \param[in] lhs The quantity on the rhs.
  //! \param[in] rhs The quantity on the rhs.
  //! \return true if equality.
  friend bool operator==(const tuple_t& lhs, const tuple_t& rhs)
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
  friend auto & operator<<(std::ostream& os, const tuple_t& a)
  {
    os << "{";
    utils::static_for_each( a.data_, 
                     [](auto & tup) { std::cout << " [ " << tup << " ]"; } );
    os << " }";
    return os;
  }
  

private:

  //! \brief The main data container, which is just a std::tuple.
  tuple_data_t data_;

};

////////////////////////////////////////////////////////////////////////////////
// Friend functions
////////////////////////////////////////////////////////////////////////////////
 
 
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





} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
