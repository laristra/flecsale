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
template <typename L, typename T, std::size_t... N> 
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

  //! the layout type
  using layout_type = L;

  //! size in each dimension 
  static constexpr size_type dimensions  = sizeof...(N);
  //! the total length of storage
  static constexpr size_type static_size = utils::multiply(N...);

private:

  //! \brief The main data container, which is just a std::array.
  T elems_[ static_size ];

  //! \brief The individual dimensions
  static constexpr size_type dims_[ dimensions ] = {N...};

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
      ( sizeof...(Args) == static_size &&
        sizeof...(Args) >= 2 ) 
    >
  >
  array(Args&&... args) : elems_{ static_cast<T>(args)... }
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
                  iterator  end()       { return elems_+static_size; }
  constexpr const_iterator  end() const { return elems_+static_size; }
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


  //! \brief return the ith element ( uses 1d index only )
  reference operator[](size_type i) 
  { 
    assert( i < static_size && "out of range" );
    return elems_[i];
  }
        
  const_reference operator[](size_type i) const 
  {     
    assert( i < static_size && "out of range" );
    return elems_[i]; 
  }

  //! \brief return the ith element ( allows multiple dimensions )
  template <typename... Args>
  std::enable_if_t< sizeof...(Args) == sizeof...(N), reference >
  operator()(Args... i) 
  { 
    assert_ranges( i... );
    auto ind = layout_type::element( i..., N... );
    return elems_[ind];
  }


  template <typename... Args>
  std::enable_if_t< sizeof...(Args) == sizeof...(N), const_reference >
  operator()(Args... i) const
  { 
    assert_ranges( i... );
    auto ind = layout_type::element( i..., N... );
    return elems_[ind];
  }


  //! \brief at() with range check
  template< typename... Args >
  std::enable_if_t< sizeof...(Args) == sizeof...(N), reference >
  at(Args... i) { 
    check_ranges(i...); 
    auto ind = layout_type::element( i..., N... );
    return elems_[ind];
  }

  template< typename... Args >
  std::enable_if_t< sizeof...(Args) == sizeof...(N), const_reference >
  at(Args... i) const 
  { 
    check_ranges(i...); 
    auto ind = layout_type::element( i..., N... );
    return elems_[ind];
  }
    
  //! \brief return the first element
  reference front() 
  { return elems_[0]; }
        
  const_reference front() const 
  { return elems_[0]; }
        
  //! \brief return the last element
  reference back() 
  { return elems_[static_size-1]; }
        
  const_reference back() const 
  {  return elems_[static_size-1]; }


  //  \brief direct access to data (read-only)
  const T* data() const { return elems_; }
  T* data() { return elems_; }

  // use array as C array (direct read/write access to data)
  T* c_array() { return elems_; }

  //===========================================================================
  // Capacity
  //===========================================================================

  //! \brief return the size
  static constexpr size_type     size() { return static_size; }
  static constexpr size_type capacity() { return size(); }

  //! \brief checks whether container is empty
  static constexpr bool empty() { return false; }
  
  //! \brief returns the maximum possible number of elements
  static constexpr size_type max_size() { return static_size; }

  //! \brief return the number of dimensions
  static constexpr size_type dims() { return dimensions; }

  //! \brief return the size in a particular dimension
  static constexpr size_type size( size_t i ) 
  { 
    assert( i<= dimensions && "dimension out of range" );
    return dims_[i];
  }


  //===========================================================================
  // operations
  //===========================================================================


  //  \brief swap contents (note: linear complexity)
  void swap (array<L,T,N...>& y) {
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
    assert( list.size() == static_size && "input list size mismatch" );
    assign( list.begin(), list.end() );
  }

  //! \brief check range (may be private because it is static)
  template< typename... Args >
  static
  std::enable_if_t< sizeof...(Args) == sizeof...(N) >
  check_ranges (Args... is) {
    utils::tuple_visit( 
                       [](auto i, auto dim) { 
                         if ( i>= dim  )
                           throw std::out_of_range("array<>: index out of range");
                       },
                       std::forward_as_tuple(is...), 
                       std::forward_as_tuple(N...) );
  }


  template< typename... Args >
  static
  std::enable_if_t< sizeof...(Args) == sizeof...(N) >
  assert_ranges ( Args... is ) 
  { 
    utils::tuple_visit( 
                       [](auto i, auto dim) { 
                         assert( i < dim && "out of range" );
                       },
                       std::forward_as_tuple(is...), 
                       std::forward_as_tuple(N...) );
  }
  
  //===========================================================================
  // Operators
  //===========================================================================


  // use std::move
  // http://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter

  //!\brief  assignment with type conversion
  template <typename T2>
  auto & operator= (const array<L,T2,N...>& rhs) {
    if ( this != &rhs )
      std::copy(rhs.begin(),rhs.end(), begin());    
    return *this;
  }

  
  //! \brief Addition binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator+=(const array<L,T2,N...> & rhs) {
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
  auto & operator-=(const array<L,T2,N...> & rhs) {
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
  auto & operator*=(const array<L,T2,N...> & rhs) {
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
  auto & operator/=(const array<L,T2,N...> & rhs) {
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
template<typename L, typename T, std::size_t... N>
bool operator==(const array<L,T,N...>& lhs, const array<L,T,N...>& rhs)
{
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template<typename L, typename T, std::size_t... N>
bool operator< (const array<L,T,N...>& x, const array<L,T,N...>& y) {
  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
}

template<typename L, typename T, std::size_t... N>
bool operator!= (const array<L,T,N...>& x, const array<L,T,N...>& y) {
  return !(x==y);
}

template<typename L, typename T, std::size_t... N>
bool operator> (const array<L,T,N...>& x, const array<L,T,N...>& y) {
  return y<x;
}
template<typename L, typename T, std::size_t... N>
bool operator<= (const array<L,T,N...>& x, const array<L,T,N...>& y) {
  return !(y<x);
}
template<typename L, typename T, std::size_t... N>
bool operator>= (const array<L,T,N...>& x, const array<L,T,N...>& y) {
  return !(x<y);
}

//! \brief  global swap(), specializes the std::swap algorithm 
template<typename L, typename T, std::size_t... N>
inline void swap (array<L,T,N...>& x, array<L,T,N...>& y) {
  x.swap(y);
}


  
//! \brief Addition operator involving two arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename L, typename T, size_t... N>
auto operator+( const array<L,T,N...>& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, typename U, size_t... N>
auto operator+( const array<L,T,N...>& lhs, 
                const U& rhs )
{
  array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e+rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator+( const U& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, size_t... N>
auto operator-( const array<L,T,N...>& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, typename U, size_t... N>
auto operator-( const array<L,T,N...>& lhs, 
                const U& rhs )
{
  array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e-rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator-( const U& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, size_t... N>
auto operator*( const array<L,T,N...>& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, typename U, size_t... N>
auto operator*( const array<L,T,N...>& lhs, 
                const U& rhs )
{
  array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e*rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator*( const U& lhs,
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, size_t... N>
auto operator/( const array<L,T,N...>& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, typename U, size_t... N>
auto operator/( const array<L,T,N...>& lhs, 
                const U& rhs )
{
  array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e/rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator/( const U& lhs, 
                const array<L,T,N...>& rhs )
{
  array<L,T,N...> tmp;
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
template <typename L, typename T, size_t... N>
auto & operator<<(std::ostream& os, const array<L,T,N...>& a)
{
  os << "(";
  for ( auto i : a ) os << " " << i;
  os << " )";
  return os;
}

//! \brief Output operator for array.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in,out] os  The ostream to dump output to.
//! \param[in]     rhs The array on the right hand side of the operator.
//! \return A reference to the current ostream.
template <typename L, typename T, std::size_t D1, std::size_t D2>
auto & operator<<(std::ostream& os, const array<L,T,D1,D2>& a)
{
  for ( std::size_t j = 0; j<D2; j++ ) { 
    os << "[";
    for ( std::size_t i = 0; i<D1; i++ ) os << " " << a(i,j);
    os << " ]" << std::endl;
  }
  return os;
}



//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename L, typename T, size_t... N>
auto dot_product(const array<L, T, N...> &a, const array<L, T, N...> &b) 
{
  return std::inner_product(a.begin(), a.end(), b.begin(), static_cast<T>(0) );
}

//! \brief Compute the magnitude of the vector
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename L, typename T, size_t... N> 
auto magnitude(const array<L, T, N...> &a) 
{
  return std::sqrt( dot_product(a,a) );
}

//! \brief Compute the magnitude of the vector
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename L, typename T, size_t... N> 
auto abs(const array<L, T, N...> &a) 
{
  return std::sqrt( dot_product(a,a) );
}

//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
template <typename L, typename T>
auto cross_product(const array<L, T, 2> &a, const array<L, T, 2> &b) 
{
  return a[0]*b[1] - a[1]*b[0];
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
