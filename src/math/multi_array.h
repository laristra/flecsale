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



struct row_major_ordering {

  using size_type = std::size_t;

  constexpr size_type stride(const size_type * ids, const size_type & N )
  {
    if ( N > 0 )
      return ids[1] * stride( &ids[1], N-1 );
    else 
      return 1;
  }
  
  template< size_type N, size_type... I >
  constexpr std::array<size_type, N> ordering_helper( const size_type (&ids)[N], std::index_sequence<I...> ) 
  {
    return {{ stride( &ids[I], N-I-1 )... }};
  }
  
  template< size_type N >
  constexpr auto operator()( const size_type (&ids)[N] )    
  {
    return ordering_helper( ids, std::make_index_sequence<N>{} );
  }
  
};

struct col_major_ordering {

  using size_type = std::size_t;

  constexpr size_type stride(const size_type * ids, const size_type & N )
  {
    if ( N > 0 )
      return ids[1] * stride( &ids[1], N-1 );
    else 
      return 1;
  }
  
  template< size_type N, size_type... I >
  constexpr std::array<size_type, N> ordering_helper( const size_type (&ids)[N], std::index_sequence<I...> ) 
  {
    return {{ stride( &ids[N-I-1], I )... }};
  }
  
  template< size_type N >
  constexpr auto operator()( const size_type (&ids)[N] )    
  {
    return ordering_helper( ids, std::make_index_sequence<N>{} );
  }
  
};

////////////////////////////////////////////////////////////////////////////////
//!  \brief The dimensioned_array type provides a general base for defining
//!  contiguous array types that have a specific dimension.
//!
//!  \tparam T The type of the array, e.g., P.O.D. type.
//!  \tparam D The dimension of the array, i.e., the number of elements
//!    to be stored in the array.
////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t... N> 
class multi_array {

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

  //! size in each dimension 
  static constexpr size_type dimensions  = sizeof...(N);
  //! the total length of storage
  static constexpr size_type elements = utils::multiply(N...);

private:

  //===========================================================================
  // Private Data
  //===========================================================================

  //! \brief The main data container, which is just a std::array.
  T elems_[ elements ];

  //! \brief The individual dimensions
  static constexpr size_type dims_[ dimensions ] = {N...};

  //! \brief The individual strides
  static constexpr std::array<size_type, sizeof...(N)> strides_ = row_major_ordering{}( {N...} );

public:

  //===========================================================================
  // Constructors / Destructors
  //===========================================================================

  //! \brief force the default constructor
  multi_array() = default;

  //! \brief force the default copy constructor
  multi_array(const multi_array &) = default;

  //!\brief fancy copy constructor with type conversion
  template <typename T2>
  multi_array(const multi_array<T2,N...>& oth) 
  {
    std::copy(oth.begin(),oth.end(), begin());    
  }

  //! \brief Constructor with one value.
  //! \param[in] val The value to set the multi_array to
  template < typename T2 >
  multi_array(const T2 & val)
  { 
    //std::cout << "multi_array (single value constructor)\n";
    fill( val ); 
  }

  //! \brief Constructor with initializer list
  //! 
  //! Initializer list is ALWAYS provided in row-major format
  //!
  //! \param[in] list the initializer list of values
  multi_array( std::initializer_list<T> list) 
  { 
    //std::cout << "multi_array (variadic constructor)\n";
    if ( list.size() == 1 ) 
      fill( *list.begin() );
    else
      assign(list);
  }
   
  //===========================================================================
  // Iterators
  //===========================================================================

  //! \brief return an iterator to the beginning of the multi_array
                  iterator  begin()       { return elems_; }
  constexpr const_iterator  begin() const { return elems_; }
  constexpr const_iterator cbegin() const { return begin(); }
        
  //! \brief return an iterator to the end of the multi_array
                  iterator  end()       { return elems_+elements; }
  constexpr const_iterator  end() const { return elems_+elements; }
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
    assert( i < elements && "out of range" );
    return elems_[i];
  }
        
  const_reference operator[](size_type i) const 
  {     
    assert( i < elements && "out of range" );
    return elems_[i]; 
  }

  //! \brief return the ith element ( allows multiple dimensions )
  template< size_type D >
  std::enable_if_t< D == dimensions, reference >
  operator[](const size_type (&ids)[D]) 
  { 
    auto ind = element( ids );
    return elems_[ind];
  }
        
  template<size_type D>
  std::enable_if_t< D == dimensions, const_reference >
  operator[](const size_type (&ids)[D]) const 
  {     
    auto ind = element( ids );
    return elems_[ind]; 
  }

  //! \brief return the ith element ( allows multiple dimensions )
  template <typename... Args>
  std::enable_if_t< sizeof...(Args) == sizeof...(N), reference >
  operator()(Args... i) 
  { 
    assert_ranges( std::forward<Args>(i)... );
    auto ind = element( std::forward<Args>(i)... );
    return elems_[ind];
  }


  template <typename... Args>
  std::enable_if_t< sizeof...(Args) == sizeof...(N), const_reference >
  operator()(Args... i) const
  { 
    assert_ranges( std::forward<Args>(i)... );
    auto ind = element( std::forward<Args>(i)... );
    return elems_[ind];
  }


  //! \brief at() with range check
  template< typename... Args >
  std::enable_if_t< sizeof...(Args) == sizeof...(N), reference >
  at(Args... i) { 
    check_ranges( std::forward<Args>(i)... ); 
    auto ind = element( std::forward<Args>(i)... );
    return elems_[ind];
  }

  template< typename... Args >
  std::enable_if_t< sizeof...(Args) == sizeof...(N), const_reference >
  at(Args... i) const 
  { 
    check_ranges( std::forward<Args>(i)... ); 
    auto ind = element( std::forward<Args>(i)... );
    return elems_[ind];
  }
    
  //! \brief return the first element
  reference front() 
  { return elems_[0]; }
        
  const_reference front() const 
  { return elems_[0]; }
        
  //! \brief return the last element
  reference back() 
  { return elems_[elements-1]; }
        
  const_reference back() const 
  {  return elems_[elements-1]; }


  //  \brief direct access to data (read-only)
  const T* data() const { return elems_; }
  T* data() { return elems_; }

  // use array as C array (direct read/write access to data)
  T* c_array() { return elems_; }

  //===========================================================================
  // Capacity
  //===========================================================================

  //! \brief return the size
  static constexpr size_type size() 
  { return elements; }

  //! \brief return the number of elements
  static constexpr size_type num_elements() 
  { return size(); }

  //! \brief return the number of dimensions
  static constexpr size_type num_dimensions() 
  { return dimensions; }

  //! \brief return the size in a particular dimension
  static constexpr const size_type * shape() 
  { return dims_; }

  //! \brief the stride associated with each array dimension
  static constexpr const size_type * strides() 
  { return strides_.data(); }


  //===========================================================================
  // operations
  //===========================================================================


  //  \brief swap contents (note: linear complexity)
  void swap (multi_array& y) {
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
    assert( list.size() == elements && "input list size mismatch" );
    assign( list.begin(), list.end() );
  }

  
  //===========================================================================
  // Operators
  //===========================================================================


  // use std::move
  // http://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter

  //!\brief  assignment with type conversion
  template <typename T2>
  auto & operator= (const multi_array<T2,N...>& rhs) {
    if ( this != &rhs )
      std::copy(rhs.begin(),rhs.end(), begin());    
    return *this;
  }

  
  //! \brief Addition binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator+=(const multi_array<T2,N...> & rhs) {
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
  auto & operator-=(const multi_array<T2,N...> & rhs) {
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
  auto & operator*=(const multi_array<T2,N...> & rhs) {
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
  auto & operator/=(const multi_array<T2,N...> & rhs) {
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
    multi_array tmp;
    std::transform( begin(), end(), tmp.begin(), std::negate<>() );    
    return tmp;
  }


  //===========================================================================
  // Utitilities
  //===========================================================================


  //! \brief unpack an array of indices and get the element index
  //! \remark this version uses an array of indices
  template< typename U, size_type D >
  static constexpr
  auto element( const U (&ids)[D] )
  {
    size_type ind = 0;
    for ( auto i=0; i<dimensions; i++ ) ind += ids[i]*strides_[i];
    return ind;
    //return std::inner_product( std::begin(ids), std::end(ids), std::begin(strides_), 0);
  }

  //! \brief compute the 1d element index
  //! \remark this version uses the variadic arguments
  template< typename... Args >
  static constexpr
  auto element( Args&&... ids )
  {
    return element( {static_cast<size_type>(ids)...} );
  }

  //! \brief check range (may be private because it is static)
  template< typename... Args >
  static
  std::enable_if_t< sizeof...(Args) == sizeof...(N) >
  check_ranges (Args... is) {
    utils::tuple_visit( 
                       [](auto i, auto dim) { 
                         if ( i>= dim  )
                           throw std::out_of_range("multi_array<>: index out of range");
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

};

////////////////////////////////////////////////////////////////////////////////
// static members
////////////////////////////////////////////////////////////////////////////////

//! \brief The individual dimensions
template <typename T, std::size_t... N> 
constexpr typename multi_array<T,N...>::size_type 
multi_array<T,N...> :: dims_[ multi_array<T,N...>::dimensions ];

//! \brief The individual dimensions
template <typename T, std::size_t... N> 
constexpr std::array<std::size_t, sizeof...(N)> multi_array<T,N...> :: strides_;

////////////////////////////////////////////////////////////////////////////////
// Friend functions
////////////////////////////////////////////////////////////////////////////////


//! \brief lexicographically compares the values in the array 
//! \param[in] lhs The quantity on the lhs.
//! \param[in] rhs The quantity on the rhs.
template<typename T, std::size_t... N>
bool operator==(const multi_array<T,N...>& lhs, const multi_array<T,N...>& rhs)
{
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template<typename T, std::size_t... N>
bool operator< (const multi_array<T,N...>& x, const multi_array<T,N...>& y) {
  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
}

template<typename T, std::size_t... N>
bool operator!= (const multi_array<T,N...>& x, const multi_array<T,N...>& y) {
  return !(x==y);
}

template<typename T, std::size_t... N>
bool operator> (const multi_array<T,N...>& x, const multi_array<T,N...>& y) {
  return y<x;
}
template<typename T, std::size_t... N>
bool operator<= (const multi_array<T,N...>& x, const multi_array<T,N...>& y) {
  return !(y<x);
}
template<typename T, std::size_t... N>
bool operator>= (const multi_array<T,N...>& x, const multi_array<T,N...>& y) {
  return !(x<y);
}

//! \brief  global swap(), specializes the std::swap algorithm 
template<typename T, std::size_t... N>
inline void swap (multi_array<T,N...>& x, multi_array<T,N...>& y) {
  x.swap(y);
}


  
//! \brief Addition operator involving two multi_arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename T, size_t... N>
auto operator+( const multi_array<T,N...>& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, typename U, size_t... N>
auto operator+( const multi_array<T,N...>& lhs, 
                const U& rhs )
{
  multi_array<T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e+rhs; } );
  return tmp;
}

template <typename T, typename U, size_t... N>
auto operator+( const U& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, size_t... N>
auto operator-( const multi_array<T,N...>& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, typename U, size_t... N>
auto operator-( const multi_array<T,N...>& lhs, 
                const U& rhs )
{
  multi_array<T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e-rhs; } );
  return tmp;
}

template <typename T, typename U, size_t... N>
auto operator-( const U& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, size_t... N>
auto operator*( const multi_array<T,N...>& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, typename U, size_t... N>
auto operator*( const multi_array<T,N...>& lhs, 
                const U& rhs )
{
  multi_array<T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e*rhs; } );
  return tmp;
}

template <typename T, typename U, size_t... N>
auto operator*( const U& lhs,
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, size_t... N>
auto operator/( const multi_array<T,N...>& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, typename U, size_t... N>
auto operator/( const multi_array<T,N...>& lhs, 
                const U& rhs )
{
  multi_array<T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e/rhs; } );
  return tmp;
}

template <typename T, typename U, size_t... N>
auto operator/( const U& lhs, 
                const multi_array<T,N...>& rhs )
{
  multi_array<T,N...> tmp;
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
template <typename T, std::size_t D1, std::size_t D2>
auto & operator<<(std::ostream& os, const multi_array<T,D1,D2>& a)
{
  for ( std::size_t j = 0; j<D2; j++ ) { 
    os << "[";
    for ( std::size_t i = 0; i<D1; i++ ) os << " " << a(i,j);
    os << " ]" << std::endl;
  }
  return os;
}


} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
