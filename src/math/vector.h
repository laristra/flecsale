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
 * \file vector.h
 * 
 * \brief Provides a dimensioned array which functions as a vector.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <array>

//! user includes
#include "ale/utils/check_types.h"

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
template <typename T, size_t D> class vector_t {

public:

  //===========================================================================
  // Typedefs
  //===========================================================================

  //! \brief the size_t type
  using size_t = std::size_t;


  //! \brief The value type.
  using value_type = T;

  //! \brief The dimension of the array.
  static constexpr size_t dimension = D;


  //===========================================================================
  // Constructors / Destructors
  //===========================================================================

  //! \brief force the default constructor
  vector_t() = default;

  //! \brief force the default copy constructor
  vector_t(const vector_t &) = default;

  //! \brief Constructor with initializer list
  //! \param[in] list the initializer list of values
  template <
    typename = typename std::enable_if< D /= 1 >::type
  >
  vector_t(std::initializer_list<T> list) {
    //std::cout << "vector_t (list constructor)\n";
    if ( list.size() == 1 )
      data_.fill(*list.begin());
    else {
      assert( list.size() == D  && " dimension size mismatch");
      std::copy(list.begin(), list.end(), data_.begin());
    }
  }

  //! \brief Constructor with initializer list
  //! \param[in] list the initializer list of values
  template <
    typename... Args,
    typename = typename std::enable_if< 
      sizeof...(Args) == D &&
      utils::are_type_t<T,Args...>::value 
    >::type
  >
  vector_t(Args&&... args)
  { 
    //std::cout << "vector_t (variadic constructor)\n";
    data_ = { args... };
  }
 
  //! \brief Constructor with one value.
  //! \param[in] val The value to set the array to
  vector_t(const T & val) 
  { 
    //std::cout << "vector_t (single value constructor)\n";
    data_.fill( val ); 
  }


  //===========================================================================
  // Accessors
  //===========================================================================


  //! \brief Return the size of the array.
  static constexpr size_t size() { return dimension; };


  //===========================================================================
  // Operators
  //===========================================================================

  //! \brief index operator
  auto & operator[](size_t i) 
  { return data_[i]; }

  //! \brief index operator (const version)
  const auto & operator[](size_t i) const 
  { return data_[i]; }


  //! \brief Assignment operator to another array.
  //! \param[in] rhs The array on the right hand side of the '='.
  //! \return A reference to the current object.
  auto & operator=(const vector_t &rhs) {
    if ( this != &rhs ) data_ = rhs.data_;
    return *this;
  }

  //! \brief Assignment operator to a constant.
  //! \param[in] val The constant on the right hand side of the '='.
  //! \return A reference to the current object.
  auto & operator=(const T &val) {
    for ( size_t i=0; i<D; i++ ) 
      data_[i] = val;    
    return *this;
  }

  // use std::move
  // http://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter
  
  //! \brief Addition binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator+=(const vector_t &rhs) {
    if ( this != &rhs ) {
      for ( size_t i=0; i<D; i++ ) 
        data_[i] += rhs[i];    
    }
    return *this;
  }

  //! \brief Addiition binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator+=(const T &val) {
    for ( size_t i=0; i<D; i++ ) 
      data_[i] += val;    
    return *this;
  }

  //! \brief Subtraction binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator-=(const vector_t &rhs) {
    if ( this != &rhs ) {
      for ( size_t i=0; i<D; i++ ) 
        data_[i] -= rhs[i];
    }
    return *this;
  }

  //! \brief Subtraction binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator-=(const T &val) {
    for ( size_t i=0; i<D; i++ ) 
      data_[i] -= val;    
    return *this;
  }


  //! \brief Multiplication binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator*=(const vector_t &rhs) {
    if ( this != &rhs ) {
      for ( size_t i=0; i<D; i++ ) 
        data_[i] *= rhs[i];    
    }
    return *this;
  }

  //! \brief Multiplication binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator*=(const T &val) {
    for ( size_t i=0; i<D; i++ ) 
      data_[i] *= val;    
    return *this;
  }

  //! \brief Division binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator/=(const vector_t &rhs) {
    if ( this != &rhs ) {
      for ( size_t i=0; i<D; i++ ) 
        data_[i] /= rhs[i];    
    }
    return *this;
  }

  //! \brief Division operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  auto & operator/=(const T &val) {
    for ( size_t i=0; i<D; i++ )
      data_[i] /= val;
    return *this;
  }

  //! \brief Equivalence operator
  //! \param[in] lhs The quantity on the rhs.
  //! \param[in] rhs The quantity on the rhs.
  //! \return true if equality.
  friend bool operator==(const vector_t& lhs, const vector_t& rhs)
  {
    if ( &lhs != &rhs ) {
      for ( size_t i=0; i<D; i++ ) 
        if ( lhs.data_[i] != rhs.data_[i] )
          return false;
    }
    return true;
  }


private:

  //! \brief The main data container, which is just a std::array.
  std::array<T, D> data_;

};


////////////////////////////////////////////////////////////////////////////////
// Friend functions
////////////////////////////////////////////////////////////////////////////////
  
//! \brief Addition operator involving two vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The vector_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator+( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp += rhs;
  return tmp;
}

//! \brief Addition operator involving one vector_t and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator+( const vector_t<T,D>& lhs, 
                const T& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp += rhs;
  return tmp;
}

template <typename T, size_t D>
auto operator+( const T& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp += rhs;
  return tmp;
}

//! \brief Subtraction operator involving two vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The vector_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator-( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

//! \brief Subtraction operator involving one vector_t and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator-( const vector_t<T,D>& lhs, 
                const T& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

template <typename T, size_t D>
auto operator-( const T& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

//! \brief Multiplication operator involving two vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The vector_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator*( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp *= rhs;
  return tmp;
}


//! \brief Multiplication operator involving one vector_t and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator*( const vector_t<T,D>& lhs, 
                const T& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

template <typename T, size_t D>
auto operator*( const T& lhs,
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

//! \brief Division operator involving two vector_ts.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The vector_t on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator/( const vector_t<T,D>& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp /= rhs;
  return tmp;
}



//! \brief Division operator involving one vector_t and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The vector_t on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t D>
auto operator/( const vector_t<T,D>& lhs, 
                const T& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp /= rhs;
  return tmp;
}

template <typename T, size_t D>
auto operator/( const T& lhs, 
                const vector_t<T,D>& rhs )
{
  vector_t<T,D> tmp(lhs);
  tmp /= rhs;
  return tmp;
}

//! \brief Output operator for vector_t.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in,out] os  The ostream to dump output to.
//! \param[in]     rhs The vector_t on the right hand side of the operator.
//! \return A reference to the current ostream.
template <typename T, size_t D>
auto & operator<<(std::ostream& os, const vector_t<T,D>& a)
{
  os << "(";
  for ( size_t i=0; i<D; i++ ) 
    os << " " << a[i];
  os << " )";
  return os;
}





//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename T, size_t D>
auto dot_product(const vector_t<T, D> &a, const vector_t<T, D> &b) {

  T sum(0);
  for (size_t d = 0; d < D; ++d)
    sum += a[d] * b[d];

  return sum;
}

template <typename T>
auto dot_product(const vector_t<T, 1> &a, const vector_t<T, 1> &b) {
  return a[0] * b[0];
}

//! \brief Compute the magnitude of the vector
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename T, size_t D> 
auto magnitude(const vector_t<T, D> &a) {
  return std::sqrt( dot_product(a,a) );
}

template <typename T> 
auto magnitude(const vector_t<T, 1> &a) {
  return std::sqrt( dot_product(a,a) );
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
