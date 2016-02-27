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


  template <typename U, size_t I, size_t... J>
  struct nested_initializer_list
  {
    using nested = typename nested_initializer_list<U, J...>::type;
    using type = std::initializer_list<nested>;
  };

  template <typename U, std::size_t I>
  struct nested_initializer_list<U, I> 
  {
    using type = std::initializer_list<U>;
  };

  template< class U, size_t... I >
  using nested_initializer_list_t = typename nested_initializer_list<U, I...>::type;



  template <typename F, size_t N>
  void multi_for( size_t (&lower_bound) [N], 
                  size_t (&upper_bound) [N],
                  F && func) {
    size_t ranges[N];
    auto numel = 1;
    for (size_t i = 0; i < N; i++) {
      ranges[i] = upper_bound[i]-lower_bound[i];
      numel *= ranges[i];
    }

    for (auto idx = 0; idx < numel; idx++) {
      //if you don't need the actual indicies, you're done

      //extract indexes
      auto idx2 = idx;
      size_t indexes[N];
      for (auto i = 0; i < ranges[i]; i++) {
        indexes[i] = idx2 % ranges[i] - lower_bound[i];
        idx2 /= ranges[i];
      }
      //do stuff
      std::forward<F>(func)(idx, indexes);
    }
  } // multi_for


////////////////////////////////////////////////////////////////////////////////
//!  \brief Define the layout of the matrix using different conventions
////////////////////////////////////////////////////////////////////////////////
struct layouts {

  //============================================================================
  //!  \brief Define the layout of the matrix using a row-major convention
  //============================================================================
  struct row_major {
    static constexpr 
    auto element( size_t i, size_t j, 
                  size_t /* size_i */, size_t size_j ) noexcept
    { return i * size_j + j; }   
  };


  //============================================================================
  //!  \brief Define the layout of the matrix using column major convention
  //============================================================================
  struct column_major {
    static constexpr 
    auto element( size_t i, size_t j, 
                  size_t size_i, size_t /* size_j */ ) noexcept
    { return j * size_i + i; }   
  };

}; // namespace layout


////////////////////////////////////////////////////////////////////////////////
//!  \brief The dimensioned_array type provides a general base for defining
//!  contiguous array types that have a specific dimension.
//!
//!  \tparam T The type of the array, e.g., P.O.D. type.
//!  \tparam D The dimension of the array, i.e., the number of elements
//!    to be stored in the array.
////////////////////////////////////////////////////////////////////////////////
template <typename L, typename T, std::size_t... N> 
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

  //! the layout type
  using layout_type = L;

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
  static constexpr size_type strides_[ dimensions ] = {N...};

  //! \brief The first element index for each dimension
  static constexpr size_type index_bases_[ dimensions ] = {N...};

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
  multi_array(const multi_array<L,T2,N...>& oth) 
  {
    std::copy(oth.begin(),oth.end(), begin());    
  }

  //! \brief Constructor with one value.
  //! \param[in] val The value to set the multi_array to
  template < typename T2 >
  multi_array(const T2 & val)
  { 
    std::cout << "multi_array (single value constructor)\n";
    fill( val ); 
  }



  class strided_iterator
  {
   public:

    strided_iterator(pointer ptr)
      : ptr_(ptr), index_(0)
    {
      for ( auto i=0; i< dimensions; i++ )
        indices_[i] = 0;
    }

    strided_iterator & operator++()
    {
      indices_[dimensions-1]++;
      for ( auto i=dimensions-1; i>0; i-- )
        if ( indices_[i] == ranges_[i] ) {
          indices_[i] = 0;
          indices_[i-1]++;
        }
      index_ = unpack_indices( indices_ );
      return *this;
    }

    strided_iterator & operator=(const strided_iterator & itr)
    {
      ptr_ = itr.ptr_;
      index_ = itr.index_;
      for ( auto i=0; i< dimensions; i++ )
        indices_[i] = itr.indices_[i];
      return *this;
    }

    reference operator*() { return ptr_[index_]; }
    pointer operator->() { return &ptr_[index_]; }

    bool operator==(const strided_iterator & itr) const
    {
      return index_ == itr.index_;
    }

    bool operator!=(const strided_iterator & itr) const
    {
      return index_ != itr.index_;
    }


   private:

    pointer ptr_;
    size_type index_;
    size_type indices_[dimensions];
    size_type ranges_[dimensions] = {N...};

  };

  //! \brief Constructor with initializer list
  //! 
  //! Initializer list is ALWAYS provided in row-major format
  //!
  //! \param[in] list the initializer list of values
  multi_array( std::initializer_list<T> list) 
  { 
    std::cout << "multi_array (variadic constructor)\n";
    if ( list.size() == 1 ) {
      fill( *list.begin() );
    }
    else if ( std::is_same_v< layout_type, layouts::row_major > ) {
      assign(list);
    }
    else {
      assert( list.size() == elements && "input list size mismatch" );
      strided_iterator it( elems_ );
      auto list_it = list.begin();
      for ( auto i=0; i<elements; i++, ++it, ++list_it)
        *it = *list_it;
    }
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
  template< std::size_t D >
  std::enable_if_t< D == dimensions, reference >
  operator[](size_type (&ids)[D]) 
  { 
    auto ind = unpack_indices( ids );
    return elems_[ind];
  }
        
  template<std::size_t D>
  std::enable_if_t< D == dimensions, const_reference >
  operator[](size_type (&ids)[D]) const 
  {     
    auto ind = unpack_indices( ids );
    return elems_[ind]; 
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
  { return dims_[0]; }

  //! \brief return the number of elements
  static constexpr size_type num_elements() 
  { return elements; }

  //! \brief return the number of dimensions
  static constexpr size_type num_dimensions() 
  { return dimensions; }

  //! \brief return the size in a particular dimension
  static constexpr const size_type * shape() 
  { return dims_; }

  //! \brief the stride associated with each array dimension
  static constexpr const size_type * strides() 
  { return strides_; }

  //! \brief the numeric index of the first element for each array dimension
  static constexpr const size_type * index_bases()
  { return index_bases_; }


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
  auto & operator= (const multi_array<L,T2,N...>& rhs) {
    if ( this != &rhs )
      std::copy(rhs.begin(),rhs.end(), begin());    
    return *this;
  }

  
  //! \brief Addition binary operator involving another array.
  //! \param[in] rhs The array on the right hand side of the operator.
  //! \return A reference to the current object.
  template <typename T2>
  auto & operator+=(const multi_array<L,T2,N...> & rhs) {
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
  auto & operator-=(const multi_array<L,T2,N...> & rhs) {
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
  auto & operator*=(const multi_array<L,T2,N...> & rhs) {
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
  auto & operator/=(const multi_array<L,T2,N...> & rhs) {
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
  // Private Utilities
  //===========================================================================

private:

  //! \brief unpack an array of indices and get the element index
  //! \remark this is the main implementation
  template< typename U, std::size_t D, std::size_t... I >
  static
  auto unpack_indices_impl( U (&ids)[D], std::index_sequence<I...> )
  {
    return layout_type::element( ids[I]..., N... );
  }

  //! \brief unpack an array of indices and get the element index
  //! \remark this is the function that should be called
  template< 
    typename U, std::size_t D, 
    typename Indices = std::make_index_sequence<D>
  >
  static
  auto unpack_indices( U (&ids)[D] )
  {
    return unpack_indices_impl( ids, Indices() );
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
// Friend functions
////////////////////////////////////////////////////////////////////////////////


//! \brief lexicographically compares the values in the array 
//! \param[in] lhs The quantity on the lhs.
//! \param[in] rhs The quantity on the rhs.
template<typename L, typename T, std::size_t... N>
bool operator==(const multi_array<L,T,N...>& lhs, const multi_array<L,T,N...>& rhs)
{
  return std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template<typename L, typename T, std::size_t... N>
bool operator< (const multi_array<L,T,N...>& x, const multi_array<L,T,N...>& y) {
  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
}

template<typename L, typename T, std::size_t... N>
bool operator!= (const multi_array<L,T,N...>& x, const multi_array<L,T,N...>& y) {
  return !(x==y);
}

template<typename L, typename T, std::size_t... N>
bool operator> (const multi_array<L,T,N...>& x, const multi_array<L,T,N...>& y) {
  return y<x;
}
template<typename L, typename T, std::size_t... N>
bool operator<= (const multi_array<L,T,N...>& x, const multi_array<L,T,N...>& y) {
  return !(y<x);
}
template<typename L, typename T, std::size_t... N>
bool operator>= (const multi_array<L,T,N...>& x, const multi_array<L,T,N...>& y) {
  return !(x<y);
}

//! \brief  global swap(), specializes the std::swap algorithm 
template<typename L, typename T, std::size_t... N>
inline void swap (multi_array<L,T,N...>& x, multi_array<L,T,N...>& y) {
  x.swap(y);
}


  
//! \brief Addition operator involving two multi_arrays.
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] lhs The array on the left hand side of the operator.
//! \param[in] rhs The array on the right hand side of the operator.
//! \return A reference to the current object.
template <typename L, typename T, size_t... N>
auto operator+( const multi_array<L,T,N...>& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator+( const multi_array<L,T,N...>& lhs, 
                const U& rhs )
{
  multi_array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e+rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator+( const U& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator-( const multi_array<L,T,N...>& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator-( const multi_array<L,T,N...>& lhs, 
                const U& rhs )
{
  multi_array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e-rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator-( const U& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator*( const multi_array<L,T,N...>& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator*( const multi_array<L,T,N...>& lhs, 
                const U& rhs )
{
  multi_array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e*rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator*( const U& lhs,
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator/( const multi_array<L,T,N...>& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto operator/( const multi_array<L,T,N...>& lhs, 
                const U& rhs )
{
  multi_array<L,T,N...> tmp;
  std::transform( lhs.begin(), lhs.end(), tmp.begin(),
                  [&rhs](auto & e) { return e/rhs; } );
  return tmp;
}

template <typename L, typename T, typename U, size_t... N>
auto operator/( const U& lhs, 
                const multi_array<L,T,N...>& rhs )
{
  multi_array<L,T,N...> tmp;
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
auto & operator<<(std::ostream& os, const multi_array<L,T,N...>& a)
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
auto & operator<<(std::ostream& os, const multi_array<L,T,D1,D2>& a)
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
auto dot_product(const multi_array<L, T, N...> &a, const multi_array<L, T, N...> &b) 
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
auto magnitude(const multi_array<L, T, N...> &a) 
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
auto abs(const multi_array<L, T, N...> &a) 
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
auto cross_product(const multi_array<L, T, 2> &a, const multi_array<L, T, 2> &b) 
{
  return a[0]*b[1] - a[1]*b[0];
}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
