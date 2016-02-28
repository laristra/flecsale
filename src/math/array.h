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
 * \file array.h
 * 
 * \brief Provides a dimensioned array which functions as a vector.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <algorithm>
#include <array>
#include <cassert>

// user includes
#include "ale/std/type_traits.h"
#include "ale/utils/check_types.h"
#include "ale/utils/template_helpers.h"
#include "ale/utils/tuple_visit.h"

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
template <typename T, std::size_t N> 
class array {

public:

  //===========================================================================
  // Typedefs
  //===========================================================================


  using value_type      = T;
  using reference       = T &;
  using pointer         = T *;
  using const_reference = const T &;
  using const_pointer   = const T *;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;

  //! iterator support
  using iterator        = pointer;
  using const_iterator  = const_pointer;
  //! reverse iterator support
  using reverse_iterator       = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  //! array size
  static constexpr size_type length  = N;

private:

  //===========================================================================
  // Private Data
  //===========================================================================

  //! \brief The main data container, which is just a std::array.
  T elems_[ length ];

public:

  //===========================================================================
  // Constructors / Destructors
  //===========================================================================

  //! \brief force the default constructor
  array() = default;

  //! \brief force the default copy constructor
  array(const array &) = default;

  //! \brief Constructor with initializer list
  //! \param[in] list the initializer list of values
  template <
    typename... Args,
    typename = std::enable_if_t< 
      ( sizeof...(Args) == N && sizeof...(Args) >= 2 )
    >
  >
  array(Args&&... args) : 
    elems_{ static_cast<T>( std::forward<Args>(args) )... }
  { 
    //std::cout << "array (variadic constructor)\n";
  }
 
  //! \brief Constructor with one value.
  //! \param[in] val The value to set the array to
  template < typename T2 >
  array(const T2 & val) 
  { 
    //std::cout << "array (single value constructor)\n";
    fill( val ); 
  }
   
  //===========================================================================
  // Iterators
  //===========================================================================

  //! \brief return an iterator to the beginning of the array
                  iterator  begin()       { return elems_; }
  constexpr const_iterator  begin() const { return elems_; }
  constexpr const_iterator cbegin() const { return begin(); }
        
  //! \brief return an iterator to the end of the array
                  iterator  end()       { return elems_+size(); }
  constexpr const_iterator  end() const { return elems_+size(); }
  constexpr const_iterator cend() const { return end(); }


  //! \brief return a reverse iterator to the beginning of the aray
  reverse_iterator rbegin() 
  { return reverse_iterator(end()); }
  
  const_reverse_iterator rbegin() const 
  { return const_reverse_iterator(end()); }
  
  const_reverse_iterator crbegin() const 
  { return const_reverse_iterator(end()); }

  //! \brief return a reverse iterator to the end of the aray
  reverse_iterator rend() 
  { return reverse_iterator(begin()); }

  const_reverse_iterator rend() const 
  { return const_reverse_iterator(begin()); }
  
  const_reverse_iterator crend() const 
  { return const_reverse_iterator(begin()); }
 
  //===========================================================================
  // Element Access
  //===========================================================================

  //! \brief return the ith element
  reference operator[](size_type i) 
  { 
    assert( i < size() && "out of range" );
    return elems_[i];
  }
        
  const_reference operator[](size_type i) const 
  {     
    assert( i < size() && "out of range" );
    return elems_[i]; 
  }

  //! \brief return the ith element
  reference operator()(size_type i) 
  { 
    assert( i < size() && "out of range" );
    return elems_[i];
  }

  const_reference operator()(size_type i) const
  { 
    assert( i < size() && "out of range" );
    return elems_[i];
  }


  //! \brief at() with range check
  reference at(size_type i) 
  { 
    return i >= size() ? 
      throw std::out_of_range("array<>: index out of range") : 
      elems_[i];
  }

  const_reference at(size_type i) const
  { 
    return i >= size() ? 
      throw std::out_of_range("array<>: index out of range") : 
      elems_[i];
  }
    
  //! \brief return the first element
  reference front() 
  { return elems_[0]; }
        
  const_reference front() const 
  { return elems_[0]; }
        
  //! \brief return the last element
  reference back() 
  { return elems_[size()-1]; }
        
  const_reference back() const 
  {  return elems_[size()-1]; }


  //  \brief direct access to data (read-only)
  const T* data() const { return elems_; }
  T* data() { return elems_; }

  // use array as C array (direct read/write access to data)
  T* c_array() { return elems_; }

  //===========================================================================
  // Capacity
  //===========================================================================

  //! \brief return the size
  static constexpr size_type     size() { return length; }
  static constexpr size_type capacity() { return size(); }

  //! \brief checks whether container is empty
  static constexpr bool empty() { return false; }
  
  //! \brief returns the maximum possible number of elements
  static constexpr size_type max_size() { return size(); }


  //===========================================================================
  // operations
  //===========================================================================

  //  \brief swap contents (note: linear complexity)
  void swap (array& y) 
  {
    std::swap(elems_, y.elems_);
  }

  //! \brief assign one value to all elements
  void fill(const T& value)
  {
    std::fill_n(begin(),size(),value);
  }

  //! \brief Replaces the contents of the container. 
  //! \tparam InputIt  The input iterator type
  //! \param [in] first  the start of the range to copy the elements from
  //! \param [in] last   the end of the range to copy the elements from
  template < class InputIt >
  void assign(InputIt first, InputIt last) 
  { 
    std::copy( first, last, begin() );
  }

  void assign( std::initializer_list<T> list ) 
  { 
    assert( list.size() == size() && "input list size mismatch" );
    assign( list.begin(), list.end() );
  }

  //===========================================================================
  // Operators
  //===========================================================================

  // use std::move
  // http://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter

  //!\brief  assignment with type conversion
  template <typename T2>
  auto & operator= (const array<T2,N>& rhs) {
    if ( this != &rhs )
      std::copy(rhs.begin(),rhs.end(), begin());    
    return *this;
  }

  
  //! \brief Addition binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator+=(const array<T2,N> & rhs) {
    std::transform( begin(), end(), rhs.begin(), 
                    begin(), std::plus<>() );
    //for ( size_type i=0; i<N; i++ ) elems_[i] += rhs.elems_[i];    
    return *this;
  }

  //! \brief Addiition binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator+=(const T2 & val) {
    std::transform( begin(), end(), begin(),
                    [&val](auto & d) { return d + val; } );
    //for ( size_type i=0; i<N; i++ ) elems_[i] += val;    
    return *this;
  }

  //! \brief Subtraction binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator-=(const array<T2,N> & rhs) {
    std::transform( begin(), end(), rhs.begin(), 
                    begin(), std::minus<>() );    
    //for ( size_type i=0; i<N; i++ ) elems_[i] -= rhs.elems_[i];    
    return *this;
  }

  //! \brief Subtraction binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator-=(const T2 & val) {
    std::transform( begin(), end(), begin(),
                    [&val](auto & d) { return d - val; } );
    //for ( size_type i=0; i<N; i++ ) elems_[i] -= val;    
    return *this;
  }


  //! \brief Multiplication binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2> 
  auto & operator*=(const array<T2,N> & rhs) {
    std::transform( begin(), end(), rhs.begin(), 
                    begin(), std::multiplies<>() );    
    //for ( size_type i=0; i<N; i++ ) elems_[i] *= rhs.elems_[i];    
    return *this;
  }

  //! \brief Multiplication binary operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator*=(const T2 & val) {
    std::transform( begin(), end(), begin(),
                    [&val](auto & d) { return d * val; } );
    //for ( size_type i=0; i<N; i++ ) elems_[i] *= val;    
    return *this;
  }

  //! \brief Division binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator/=(const array<T2,N> & rhs) {
    std::transform( begin(), end(), rhs.begin(), 
                    begin(), std::divides<>() );    
    //for ( size_type i=0; i<N; i++ ) elems_[i] /= rhs.elems_[i];    
    return *this;
  }

  //! \brief Division operator involving a constant.
  //! \param[in] val The constant on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator/=(const T2 & val) {
    std::transform( begin(), end(), begin(),
                    [&val](auto & d) { return d / val; } );
    //for ( size_type i=0; i<N; i++ ) elems_[i] /= val;
    return *this;
  }

  //! \brief Unary - operator.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  auto operator-() const {
    array tmp;
    std::transform( begin(), end(), tmp.begin(), std::negate<>() );    
    return tmp;
  }

};


////////////////////////////////////////////////////////////////////////////////
// Friend functions
////////////////////////////////////////////////////////////////////////////////


//! \brief lexicographically compares the values in the array 
//! \param[in] lhs The quantity on the lhs.
//! \param[in] rhs The quantity on the rhs.
template<typename T, std::size_t N>
bool operator==(const array<T,N>& lhs, const array<T,N>& rhs)
{
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template<typename T, std::size_t N>
bool operator< (const array<T,N>& x, const array<T,N>& y) {
  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
}

template<typename T, std::size_t N>
bool operator!= (const array<T,N>& x, const array<T,N>& y) {
  return !(x==y);
}

template<typename T, std::size_t N>
bool operator> (const array<T,N>& x, const array<T,N>& y) {
  return y<x;
}
template<typename T, std::size_t N>
bool operator<= (const array<T,N>& x, const array<T,N>& y) {
  return !(y<x);
}
template<typename T, std::size_t N>
bool operator>= (const array<T,N>& x, const array<T,N>& y) {
  return !(x<y);
}

//! \brief  global swap(), specializes the std::swap algorithm 
template<typename T, std::size_t N>
inline void swap (array<T,N>& x, array<T,N>& y) {
  x.swap(y);
}


  
//! \brief Addition operator involving two arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, std::size_t N>
auto operator+( const array<T,N>& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), 
                  tmp.begin(), std::plus<>() );    
  return tmp;
}

//! \brief Addition operator involving one array and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, typename U, std::size_t N>
auto operator+( const array<T,N>& lhs, 
                const U& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e+rhs; } );
  return tmp;
}

template <typename T, typename U, std::size_t N>
auto operator+( const U& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(),
                  [&lhs](auto & e) { return lhs+e; } );
  return tmp;
}

//! \brief Subtraction operator involving two arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, std::size_t N>
auto operator-( const array<T,N>& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), 
                  tmp.begin(), std::minus<>() );    
  return tmp;
}

//! \brief Subtraction operator involving one array and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, typename U, std::size_t N>
auto operator-( const array<T,N>& lhs, 
                const U& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e-rhs; } );
  return tmp;
}

template <typename T, typename U, std::size_t N>
auto operator-( const U& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(),
                  [&lhs](auto & e) { return lhs-e; } );
  return tmp;
}

//! \brief Multiplication operator involving two arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, std::size_t N>
auto operator*( const array<T,N>& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), 
                  tmp.begin(), std::multiplies<>() );    
  return tmp;
}


//! \brief Multiplication operator involving one array and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, typename U, std::size_t N>
auto operator*( const array<T,N>& lhs, 
                const U& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e*rhs; } );
  return tmp;
}

template <typename T, typename U, std::size_t N>
auto operator*( const U& lhs,
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(),
                  [&lhs](auto & e) { return lhs*e; } );
  return tmp;
}

//! \brief Division operator involving two arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, std::size_t N>
auto operator/( const array<T,N>& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), 
                  tmp.begin(), std::divides<>() );    
  return tmp;
}



//! \brief Division operator involving one array and a scalar.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The scalar on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, typename U, std::size_t N>
auto operator/( const array<T,N>& lhs, 
                const U& rhs )
{
  array<T,N> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e/rhs; } );
  return tmp;
}

template <typename T, typename U, std::size_t N>
auto operator/( const U& lhs, 
                const array<T,N>& rhs )
{
  array<T,N> tmp;
  std::transform( rhs.begin(), rhs.end(), tmp.begin(),
                  [&lhs](auto & e) { return lhs/e; } );
  return tmp;
}

//! \brief Output operator for array.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in,out] os  The ostream to dump output to.
//! \param[in]     rhs The array on the right hand side of the operator.
//! \return A reference to the current ostream.
template <typename T, std::size_t N>
auto & operator<<(std::ostream& os, const array<T,N>& a)
{
  os << "(";
  for ( auto i : a ) os << " " << i;
  os << " )";
  return os;
}

//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename T>
auto cross_product(const array<T, 2> &a, const array<T, 2> &b) 
{
  return a[0]*b[1] - a[1]*b[0];
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
