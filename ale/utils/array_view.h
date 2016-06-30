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
 * \file array_ref.h
 * 
 * \brief A reference array to avoid a million overloads.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/std/type_traits.h"
#include "ale/utils/type_traits.h"

namespace ale {
namespace utils {

namespace detail {

////////////////////////////////////////////////////////////////////////////////    
/// \brief A Helper to identify if this is a container
////////////////////////////////////////////////////////////////////////////////    
template<typename T, typename _ = void>
struct is_container : std::false_type {};

template<typename... Ts>
struct is_container_helper {};

template<typename T>
struct is_container<
  T,
  std::conditional_t<
    false,
    is_container_helper<
      decltype(std::declval<T>().size()),
      decltype(std::declval<T>().data())
    >,
    void
  >
> : public std::true_type {};

template< typename T >
constexpr bool is_container_v = is_container<T>::value;


} // namespace detail


////////////////////////////////////////////////////////////////////////////////    
//! forward declarations
////////////////////////////////////////////////////////////////////////////////    
template<int Rank>
class bounds_iterator;


////////////////////////////////////////////////////////////////////////////////    
/// \brief the index represenation
/// An N-dimensional vector identifying a point in N-dimensional space
////////////////////////////////////////////////////////////////////////////////    
template <int Rank>
class index_t {
public:

  //============================================================================
  /// \name types
  //============================================================================
  /// @{
  //! \brief the number of dimensions in the rectangle
  static constexpr auto rank = Rank;
  //! \brief the type of index, which is ptrdiff_t
  //! This allows the index to both address every byte in the largest allocation 
  //! and express any offset in the same
  using value_type = std::ptrdiff_t;
  //! \brief the reference to a value
  using reference = value_type &;
  //! \brief the constant reference to a value
  using const_reference = const value_type &;
  //! \brief the reference to a value
  using pointer = value_type *;
  //! \brief the constant reference to a value
  using const_pointer = const value_type *;
  //! \brief the size type used 
  using size_type = std::size_t;
  //! \brief the iterator types
  using iterator = pointer;
  //! \brief the constant interator
  using const_iterator = const_pointer;
  /// @}

  //============================================================================
  /// \name construct/copy
  //============================================================================
  /// @{
  
  //! \brief the default constructor
  constexpr index_t() noexcept = default;
  //! \brief the default copy constructor
  constexpr index_t(const index_t &) noexcept = default;
  //! \brief the default move constructor
  constexpr index_t(index_t &&) noexcept = default;
  //! \brief the default assignment operator
  index_t& operator=(const index_t &) noexcept = default;
  //! \brief the default move assignment operator
  index_t& operator=(index_t &&) noexcept = default;


  //! \brief constructor with values, must be the same number as rank
  template <
    typename... T,
    typename = typename std::enable_if_t< (sizeof...(T) == Rank) >
  > 
  constexpr index_t(T... ts) noexcept : 
    index_{ static_cast<value_type>( std::forward<T>(ts) )...} 
  {};


  //! \brief constructor with static array, must be the same size as rank\
  //! \tparam N the rank of the array, must equal Rank
  //! \param [in] vals  the values to set
  template <
    int N,
    typename = typename std::enable_if_t< (N == Rank) >
  > 
  constexpr index_t(const value_type (&vals)[N]) noexcept
  {
    for ( size_type i=0; i<rank; i++) index_[i] = vals[i];
  };
  /// @}
  

  //============================================================================
  /// \name Access components
  //============================================================================
  /// @{
  //! \brief access the individual bounds for a specific dimension
  //! \param [in] component_idx the component to return
  //! \return the bounds
  reference operator[](size_type component_idx) noexcept 
  { return index_[component_idx]; };

  //! \brief access the individual bounds for a specific dimension
  //! \param [in] component_idx the component to return
  //! \return the bounds
  constexpr const_reference operator[](size_type component_idx) const noexcept 
  { return index_[component_idx]; };
  /// @}

  //! \brief return the beginning iterator
  iterator begin() noexcept 
  { return &index_[0]; }
  const_iterator begin() const noexcept 
  { return &index_[0]; }

  //! \brief return the end iterator
  iterator end() noexcept 
  { return &index_[rank]; }
  const_iterator end() const noexcept 
  { return &index_[rank]; }

  //============================================================================
  /// \name Comparison operators
  //============================================================================
  /// @{
  //! \brief test for equivalency
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return true if equal
  bool operator==(const index_t & rhs) const noexcept 
  {  
    return std::equal( 
      std::begin(index_), std::end(index_), std::begin(rhs.index_) 
    ); 
  };

  //! \brief test for non-equivalency
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return true if not equal
  bool operator!=(const index_t & rhs) const noexcept
  {
    return !std::equal( 
      std::begin(index_), std::end(index_), std::begin(rhs.index_) 
    ); 
  };
  /// @}

  //============================================================================
  /// \name Arithmatic operators
  //============================================================================
  /// @{

  //! \brief addition operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return return the result of the operation
  index_t operator+(const index_t<rank>& rhs) const noexcept
  { 
    index_t tmp;
    for ( size_type i=0; i<rank; i++ )
      tmp[i] = this->operator[](i) + rhs[i];
    return tmp;
  };

  //! \brief minus operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return return the result of the operation
  index_t operator-(const index_t<rank>& rhs) const noexcept
  { 
    index_t tmp;
    for ( size_type i=0; i<rank; i++ )
      tmp[i] = this->operator[](i) - rhs[i];
    return tmp;
  };

  //! \brief multiply operator for arithmetic types only
  //! \tparam T must be an arithmetic type
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return return the result of the operation
  template< typename T >
  std::enable_if_t<std::is_arithmetic<T>::value, index_t> 
  operator*(T v) const noexcept
  { 
    index_t tmp;
    for ( size_type i=0; i<rank; i++ )
      tmp[i] = this->operator[](i) * v;
    return tmp;
  };
  
  //! \brief division operator for arithmetic types only
  //! \tparam T must be an arithmetic type
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return return the result of the operation
  template< typename T >
  std::enable_if_t<std::is_arithmetic<T>::value, index_t> 
  operator/(T v) const noexcept
  { 
    index_t tmp;
    for ( size_type i=0; i<rank; i++ )
      tmp[i] = this->operator[](i) / v;
    return tmp;
  };

  //! \brief in-place addition operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return a reference to this object
  index_t& operator+=(const index_t<rank>& rhs) noexcept
  { 
    for ( size_type i=0; i<rank; i++ ) this->operator[](i) += rhs[i];
    return *this;
  };

  //! \brief in-place minus operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return a reference to this object
  index_t& operator-=(const index_t<rank>& rhs) noexcept
  { 
    for ( size_type i=0; i<rank; i++ ) this->operator[](i) -= rhs[i];
    return *this;
  };

  //! \brief in-place multiply operator for arithmetic types only
  //! \tparam T must be an arithmetic type
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return a reference to this object
  template< typename T >
  std::enable_if_t<std::is_arithmetic<T>::value, index_t&> 
  operator*=(T v) noexcept
  { 
    for ( size_type i=0; i<rank; i++ ) this->operator[](i) *= v;
    return *this;
  };

  //! \brief in-place division operator for arithmetic types only
  //! \tparam T must be an arithmetic type
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return a reference to this object
  template< typename T >
  std::enable_if_t<std::is_arithmetic<T>::value, index_t&> 
  operator/=(T v) noexcept
  { 
    for ( size_type i=0; i<rank; i++ ) this->operator[](i) /= v;
    return *this;
  };

  /// @}

  //============================================================================
  /// \name Unary operators
  //============================================================================
  /// @{

  //! \brief unary + operator
  //! \return a copy of this object
  auto operator+() const noexcept
  { return *this; }

  //! \brief unary - operator
  //! \return a copy of this object negated
  auto operator-() const noexcept
  {
    index_t tmp(*this);
    for ( size_type i=0; i<rank; i++ ) tmp[i] = -tmp[i];
    return tmp;
  }

  //! \brief prefix increement
  //! \remark only valid for Rank=1
  //! \tparam R the rank
  //! \return a reference to the incremented object
  template< int R = rank >
  std::enable_if_t< (R==1), index_t& > operator++() noexcept
  { 
    ++index_[0]; 
    return *this;
  }

  //! \brief prefix decrement
  //! \remark only valid for Rank=1
  //! \tparam R the rank
  //! \return a reference to the decremented object
  template< int R = rank >
  std::enable_if_t< (R==1), index_t& > operator--() noexcept
  { 
    --index_[0]; 
    return *this;
  }

  //! \brief postfix increement
  //! \remark only valid for Rank=1
  //! \tparam R the rank
  //! \return a copy of the incremented object
  template< int R = rank >
  std::enable_if_t< (R==1), index_t > operator++(int) noexcept
  { 
    auto tmp = *this;
    ++(*this); 
    return tmp;
  }

  //! \brief postfix decrement
  //! \remark only valid for Rank=1
  //! \tparam R the rank
  //! \return a copy of the decremented object
  template< int R = rank >
  std::enable_if_t< (R==1), index_t > operator--(int) noexcept
  { 
    auto tmp = *this;
    --(*this); 
    return tmp;
  }
  /// @}

  //============================================================================
  /// \name Friend operators
  //============================================================================
  /// @{


  //! \brief namespace addition operator
  //! \tparam R the rank of the bounds object
  //! \param [in] lhs the object on the left hand side of the operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return the result of the operation
  template< int R >
  friend index_t<R> operator+(const index_t<R>& lhs, const index_t<R>& rhs) noexcept;

  //! \brief namespace minus operator
  //! \tparam R the rank of the bounds object
  //! \param [in] lhs the object on the left hand side of the operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return the result of the operation
  template< int R >
  friend index_t<R> operator-(const index_t<R>& lhs, const index_t<R>& rhs) noexcept;

  //! \brief namespace multiplication operator
  //! \tparam R the rank of the bounds object
  //! \param [in] lhs the object on the left hand side of the operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return the result of the operation
  template< int R, typename T > 
  friend index_t<R> operator*( T lhs, const index_t<R>& rhs ) noexcept;

  //! \brief namespace division operator
  //! \tparam R the rank of the bounds object
  //! \param [in] lhs the object on the left hand side of the operator
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return the result of the operation
  template< int R, typename T > 
  friend index_t<R> operator/( const index_t<R>& lhs, T rhs ) noexcept;

  /// @}

protected:

  //============================================================================
  /// \name private data
  //============================================================================
  /// @{  
  value_type index_[rank] = {};
  /// @}

};


//==============================================================================
/// \brief friend operators
//==============================================================================
/// @{

template< int Rank >
index_t<Rank> operator+(const index_t<Rank>& lhs, const index_t<Rank>& rhs) noexcept
{
  index_t<Rank> tmp(lhs);
  tmp += rhs;
  return tmp;
}

template< int Rank >
index_t<Rank> operator-(const index_t<Rank>& lhs, const index_t<Rank>& rhs) noexcept
{
  index_t<Rank> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

template< int Rank, typename T > 
std::enable_if_t< std::is_arithmetic<T>::value, index_t<Rank> >
operator*( T lhs, const index_t<Rank>& rhs ) noexcept
{
  index_t<Rank> tmp(rhs);
  tmp *= lhs;
  return tmp;
}

template< int Rank, typename T > 
std::enable_if_t< std::is_arithmetic<T>::value, index_t<Rank> >
operator/( const index_t<Rank>& lhs, T rhs ) noexcept
{
  index_t<Rank> tmp(lhs);
  tmp /= rhs;
  return tmp;
}

/// @}

////////////////////////////////////////////////////////////////////////////////    
/// \brief the bounds represenation
///
/// This is an N-dimensional axis aligned rectangle with the minimum point at 0
////////////////////////////////////////////////////////////////////////////////    
template <int Rank>
class bounds_t {

public:

  //============================================================================
  /// \name types
  //============================================================================
  /// @{
  //! \brief the number of dimensions in the rectangle
  static constexpr auto rank = Rank;
  //! \brief the index type
  using index_type = index_t<Rank>;
  //! \brief the type of index, which is ptrdiff_t
  //! This allows the index to both address every byte in the largest allocation 
  //! and express any offset in the same
  using value_type = typename index_type::value_type;
  //! \brief the reference to a value
  using reference = value_type &;
  //! \brief the constant reference to a value
  using const_reference = const value_type &;
  //! \brief the size type used 
  using size_type = std::size_t;
  /// @}

  //============================================================================
  /// \name construct/copy
  //============================================================================
  /// @{
  
  //! \brief the default constructor
  constexpr bounds_t() noexcept = default;
  //! \brief the default copy constructor
  constexpr bounds_t(const bounds_t &) noexcept = default;
  //! \brief the default move constructor
  constexpr bounds_t(bounds_t &&) noexcept = default;
  //! \brief the default assignment operator
  bounds_t& operator=(const bounds_t &) noexcept = default;
  //! \brief the default move assignment operator
  bounds_t& operator=(bounds_t &&) noexcept = default;

  //! \brief constructor with values, must be the same number as rank
  template <
    typename... T,
    typename = typename std::enable_if_t< (sizeof...(T) == Rank) >
  > 
  constexpr bounds_t(T... ts) noexcept : 
    bounds_{ std::forward<T>(ts)...} 
  {
    set_strides();
  };


  //! \brief constructor with static array, must be the same size as rank\
  //! \tparam N the rank of the array, must equal Rank
  //! \param [in] vals  the values to set
  template <
    int N,
    typename = typename std::enable_if_t< (N == Rank) >
  > 
  constexpr bounds_t(const value_type (&vals)[N]) noexcept :
    bounds_(vals)
  {
    set_strides();
  };
    

  /// @}

  //============================================================================
  /// \name Access components
  //============================================================================
  /// @{
  //! \brief access the individual bounds for a specific dimension
  //! \param [in] component_idx the component to return
  //! \return the bounds
  reference operator[](size_type component_idx) noexcept 
  { return bounds_[component_idx]; };

  //! \brief access the individual bounds for a specific dimension
  //! \param [in] component_idx the component to return
  //! \return the bounds
  constexpr const_reference operator[](size_type component_idx) const noexcept 
  { return bounds_[component_idx]; };
  /// @}


  //============================================================================
  /// \name Comparison operators
  //============================================================================
  /// @{
  //! \brief test for equivalency
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return true if equal
  bool operator==(const bounds_t & rhs) const noexcept 
  {  
    return std::equal( 
      bounds_.begin(), bounds_.end(), rhs.bounds_.begin() 
    ); 
  };

  //! \brief test for non-equivalency
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return true if not equal
  bool operator!=(const bounds_t & rhs) const noexcept
  {
    return !std::equal( 
      bounds_.begin(), bounds_.end(), rhs.bounds_.end() 
    ); 
  };
  /// @}

  //============================================================================
  /// \name bounds_t functions
  //============================================================================
  /// @{

  //! \brief the volume of the domain
  //! \return the product of all dimensions
  constexpr auto size() const noexcept
  { 
    return std::accumulate(
      bounds_.begin(), bounds_.end(), 
      static_cast<size_type>(1), std::multiplies<size_type>()
    );    
  }

  //! \brief is particular index inside the domain bounds
  //! \param [in] idx the index to check against the bounds
  //! \return true if idx is within the bounds
  bool contains(const index_t<rank>& idx) const noexcept
  {
    for ( size_type i=0; i<rank; i++ )
      if ( idx[i] < 0 || idx[i] >= this->operator[](i) )
        return false;
    return true;
  }
  bool contains(const bounds_t& bnds) const noexcept
  {
    for ( size_type i=0; i<rank; i++ )
      if ( bnds[i] < 0 || bnds[i] > this->operator[](i) )
        return false;
    return true;
  }

  //! \brief return the stride of the array
  constexpr const index_type & stride() const noexcept
  { return stride_; }

  //! \brief comnpute a 1d index from a multidimensional one
  //! \param [in] idx the index to linearize
  constexpr auto linearize( const index_type & idx ) const noexcept {
    auto offset = idx[0] * stride_[0];
    for ( size_type i=1; i<rank; i++ ) offset += idx[i] * stride_[i];
    return offset; 
  }

  //! \brief return the begining iterator {0} 
  //! \return an iterator to the first index
  auto begin() const noexcept 
  { return bounds_iterator<rank>( *this ); };

  //! \brief return the end iterator which is just past the bounds
  //! \return an iterator past the bounds
  auto end() const noexcept
  { return bounds_iterator<rank>( *this, index_t<rank>(this->bounds_) ); }
  /// @}

  //! \brief set the stride using the default scheme
  constexpr void set_strides() noexcept
  { 
    set_row_major();
  };  

protected:

  //============================================================================
  /// \name helper functions
  //============================================================================

  //! \brief set the stride to row major based on the bounds
  constexpr void set_row_major() noexcept
  { 
    stride_[rank-1] = 1;
    for ( int i=rank-2; i>-1; i-- )
      stride_[i] = stride_[i+1]*bounds_[i+1];
  };  

  //! \brief set the stride to column major based on the bounds
  constexpr void set_column_major() noexcept
  { 
    stride_[0] = 1;
    for ( int i=1; i<rank; i++ )
      stride_[i] = stride_[i-1]*bounds_[i-1];
  };  


  //============================================================================
  /// \name private data
  //============================================================================
  /// @{  
  index_type bounds_;
  //! \brief the strides of the data
  index_type stride_;

  /// @}
  
};

////////////////////////////////////////////////////////////////////////////////    
/// \brief the bounds represenation
////////////////////////////////////////////////////////////////////////////////    
template <int Rank>
class bounds_iterator {
public:

  //============================================================================
  /// \name types
  //============================================================================
  /// @{

  //! \brief an alias for the bounds type
  using bounds_type = bounds_t<Rank>;
  //! \brief an alias for the index type
  using index_type = typename bounds_type::index_type;

  //! \brief this is a random access iterator
  using iterator_category = std::random_access_iterator_tag;
  //! \brief the value type is an index
  using value_type = const index_type;
  //! \brief the standard different type
  using difference_type = std::ptrdiff_t;
  //! \brief a pointer to the index
  using pointer = const index_type*;
  //! \brief a reference to the index type
  using reference = const index_type;

  /// @}

  //============================================================================
  /// \name construct/copy
  //============================================================================
  /// @{

  //! \brief an explcit constructor to create a new iterator
  //! \param [in] bnd the bounds the iterator covers
  //! \param [in] curr the current value the iterator points to
  explicit bounds_iterator(
    const bounds_type & bnd, const index_type & curr = index_type{} ) noexcept : 
    bounds_(bnd), current_(curr) 
  {};

  //! \brief default copy constructor
  bounds_iterator(const bounds_iterator& rhs) = default;
  //! \brief default move constructor
  bounds_iterator(bounds_iterator&& rhs) = default;
  //! \brief default assignement operator
  bounds_iterator& operator=(const bounds_iterator& rhs) = default;
  /// @}

  //============================================================================
  /// \name Comparison operators
  //============================================================================
  /// @{
  //! \brief test for equivalency
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return true if equal
  bool operator==(const bounds_iterator & rhs) const noexcept 
  {  
    return (current_ == rhs.current_);
  };

  //! \brief test for non-equivalency
  //! \param [in] rhs the object on the right hand side of the operator
  //! \return true if not equal
  bool operator!=(const bounds_iterator & rhs) const noexcept
  {
    return (current_ != rhs.current_);
  };
  /// @}

private:

  //============================================================================
  /// \name private data
  //============================================================================
  /// @{

  //! \brief storage for the bounds
  bounds_type bounds_;
  //! \brief storage for the current index
  index_type current_;

  /// @}

};


////////////////////////////////////////////////////////////////////////////////    
/// \brief represent multidimensional views onto regular collections of 
/// objects of a uniform type.
/// 
/// The view semantics convey that objects of these types do not store the 
/// actual data, but instead enables patterns congruent to that of random 
/// access iterators or pointers. 
////////////////////////////////////////////////////////////////////////////////    
template < typename T, int Rank = 1 >
class array_view {

public:

  //============================================================================
  /// \name types
  //============================================================================

  /// @{
  //! \brief the rank of the array view
  static constexpr auto rank = Rank;
  //! \brief the index type
  using index_type = index_t<rank>;
  //! \brief the bounds type
  using bounds_type = bounds_t<rank>;
  //! \brief the size type
  using size_type = typename bounds_type::size_type;

  //! \brief the value type
  using value_type = T;
  //! \brief a pointer to the data type
  using pointer = typename std::add_pointer_t<value_type>;
  //! \brief a reference to the data type
  using reference = typename std::add_lvalue_reference_t<value_type>;
  /// @}
        
  //============================================================================
  /// \name construct/copy
  //============================================================================
  /// @{

  //! \brief default constructor
  constexpr array_view() noexcept = default;
  //! \brief default copy constructor
  constexpr array_view(const array_view&) noexcept = default; 
  //! \brief default move constructor
  constexpr array_view(array_view&&) noexcept = default;
  //! \brief default assignement operator
  array_view& operator=(const array_view&) noexcept = default;
  //! \brief default move assignement operator
  array_view& operator=(array_view&&) noexcept = default;


  //! \brief single parameter array_view constructor
  //! \remark Data must be convertible to pointer
  //! \remark This serves a copy constructor but also as an implicit conversion
  //!         between related array_view types, .e.g. converting non-const 
  //!         to const data
  //! \param [in] rhs  the array_vew to set the new view to
  template<typename ViewValueType>
  constexpr array_view( 
    const array_view<ViewValueType, rank>& rhs,
    std::enable_if_t< std::is_convertible< std::add_pointer_t<ViewValueType>, pointer>::value >* = nullptr
  ) noexcept :
    data_( rhs.data() ), bounds_( rhs.bounds() )
  {
    std::cout << "copy constructor" << std::endl;
  };

  //! \brief array_view copy assignemnt operator
  //! \remark Data must be convertible to pointer
  //! \remark This serves a copy constructor but also as an implicit conversion
  //!         between related array_view types, .e.g. converting non-const 
  //!         to const data
  //! \param [in] rhs  the array_vew to set the new view to
  template<typename ViewValueType>
  std::enable_if_t< 
    std::is_convertible< std::add_pointer_t<ViewValueType>, pointer>::value, 
    array_view&  
  >
  operator=(const array_view<ViewValueType, rank>& rhs) noexcept
  {
    data_ = rhs.data();
    bounds_ = rhs.bounds();
    std::cout << "copy assignement" << std::endl;
  }

  //! \brief single parameter container
  //! \remark only allowed for rank = 1
  //! \remark Note this constructor also allows for converting array_views with 
  //!         rank>1 to array_views with rank=  1, i.e. flattening
  //! \remark only for stl-like containers
  //! \param [in] cont  the container to set the view to
  template<
    typename Container,
    typename = typename std::enable_if_t< 
      detail::is_container_v<Container> && (rank == 1)
    >
  >
  constexpr explicit array_view(Container & cont) noexcept :
    data_( cont.data() ), 
    bounds_( static_cast<typename bounds_type::value_type>( cont.size() ) )
  {
    std::cout << "container constructor" << std::endl;
  }


  //! \brief two parameter constructor with bounds and container
  //! \remark 
  //! \remark only for stl-like containers
  //! \param [in] bounds  the bounds of the view to set
  //! \param [in] cont  the container to set the view to
  template<
    typename Container,
    typename = typename std::enable_if_t< 
      detail::is_container_v<Container>
    >
  >
  constexpr array_view(bounds_type bounds, Container& cont) noexcept :
    data_( cont.data() ), 
    bounds_( bounds )
  {
    // the container size() must be greater than or equal to the bounds size
    assert( cont.size() >= bounds.size() && "container size is smaller than bounds" );
    std::cout << "container with bounds constructor" << std::endl;
  }

  //! \brief single parameter array constructor
  //! \remark The array bounds are derived from the extents
  //! \remark only for arrays that can be converted to pointer with rank<> == rank
  //! \param [in] arr  the container to set the view to
  template<
    typename ArrayType,
    typename std::enable_if_t< 
      // can convert pointers
      std::is_convertible< std::add_pointer_t< std::remove_all_extents_t<ArrayType> >, pointer >::value
      // has same rank
      && std::rank<ArrayType>::value == rank 
      // has same type
      &&  std::is_same_v< std::remove_all_extents_t<ArrayType>, value_type >
    >* = nullptr
  >
  constexpr explicit array_view(ArrayType & arr) noexcept :
    array_view( arr, std::make_index_sequence< std::rank<ArrayType>::value >{} )
  {
    std::cout << "array constructor" << std::endl;
  }


  //! \brief two parameter constructor with bounds and pointer
  //! \remark  pointed to storage contains at least as many adjacent objects 
  //!          as the bounds size()
  //! \remark only for arrays that can be converted to pointer with rank<> == rank
  //! \param [in] bounds  the bounds of the view to set
  //! \param [in] data  the data to set the view to
  constexpr array_view(bounds_type bounds, pointer data) noexcept :
    data_( data ), bounds_( bounds )
  {
    std::cout << "pointer constructor with bounds" << std::endl;
  }
        
  /// @}
        
        
  //============================================================================
  /// \name observers
  //============================================================================
  /// @{

  //! \brief return the total size of the array
  constexpr size_type size() const { return bounds().size(); }
  constexpr size_type max_size() const {
    return std::numeric_limits<size_type>::max() / sizeof(T);
  }

  //! \brief return true if empty
  constexpr bool empty() const { return size() == 0; }

  //! \brief return the bounds of the array
  constexpr const bounds_type & bounds() const noexcept
  { return bounds_; }

  //! \brief return the stride of the array
  constexpr const index_type & stride() const noexcept
  { return bounds_.stride(); }
  /// @}
        
  //============================================================================
  /// \name element access
  //============================================================================
  /// @{

  //! \brief access the element at a specified index
  //! \remark this version has no bounds checking
  //! \param [in] idx the index
  //! \return a reference to the indexed location
  constexpr reference operator[](const index_type &idx) const noexcept
  { 
    return data_[ bounds_.linearize(idx) ]; 
  }

  //! \brief access the element at a specified index
  //! \remark this version has bounds checking
  //! \param [in] idx the index
  //! \return a reference to the indexed location
  constexpr reference at(const index_type &idx) const noexcept 
  {
    // This makes at() constexpr as long as the argument is within the
    // bounds of the array_view.    
    return bounds().contains(idx) ? this->operator[](idx) : 
      throw std::out_of_range("at() argument out of range");
  }

  template <typename FirstIndex, typename... OtherIndices>
  constexpr reference operator()(FirstIndex index, OtherIndices... indices)
  {
    index_type idx{index, indices...};
    return this->operator[](idx);
  }        
  /// \returns A pointer such that [<code>data()</code>,<code>data() +
  /// size()</code>) is a valid range. For a non-empty array_view,
  /// <code>data() == &front()</code>.
  constexpr pointer data() const noexcept { return data_; }
  /// @}

  //============================================================================
  /// \name slicing
  //============================================================================
  /// @{
                                                                 
  //! \brief slice an array, i.e. fix the most significant dimension
  //! \remark this version has no bounds checking
  //! \param [in] slice the index
  //! \return a reference to the indexed location
  template< int R = rank-1 >
  std::enable_if_t< (R>0), array_view< value_type, R > > 
  operator[]( typename index_type::value_type slice ) const noexcept
  { 
    return array_view< value_type, R >( bounds_.slice(), data_ );
  }

  //! \brief slice an array, i.e. fix the most significant dimension
  //! \remark this version has bounds checking
  //! \param [in] slice the index
  //! \return a reference to the indexed location
  template< int R = rank-1 >
  std::enable_if_t< (R>0), array_view< value_type, R > > 
  at( typename index_type::value_type slice ) const noexcept
  { 
    // This makes at() constexpr as long as the argument is within the
    // bounds of the array_view.    
    return (slice < 0 || slice >= bounds()[0] ) ? 
      throw std::out_of_range("at() argument out of range") :
      this->operator[](slice);
  }

                                                                 
  //! \brief section an array, i.e. create a sub-view
  //! \param [in] origin   the index, if non specified, use {0} origin
  //! \param [in] bnds     the new bounds, if non is specified, keep old bounds
  //! \return a the new view
  array_view section() const noexcept
  { return *this; }

  array_view section( const index_type & origin ) const noexcept
  { 
    assert( bounds().contains(origin) );
    auto bnds = bounds() - origin;
    return array_view( bnds, &this->operator[](origin) );
  }

  array_view section( 
    const index_type & origin, 
    const bounds_type& bnds ) const noexcept
  { 
    assert( bnds.contains(origin) );
    assert( bounds().contains(bnds) );
    return array_view( bnds, &this->operator[](origin) );
  }

  
  /// @}

  //============================================================================
  /// \name mutators
  //============================================================================
  /// @{
        
  /// \par Effects:
  /// Resets *this to its default-constructed state.
  void clear() { *this = array_view(); }


  //! \brief declare other variations of array_view as a friend
  template< typename U,  int R >
  friend class array_view;

protected:


  //============================================================================
  /// \name private helpers
  //============================================================================

  //! \brief constructor with static array, must be the same size as rank
  //! \tparam N the rank of the array, must equal Rank
  //! \param [in] vals  the values to set
  template < typename ArrayType, std::size_t... I > 
  constexpr array_view( ArrayType & arr, std::index_sequence<I...> ) noexcept : 
    data_( reinterpret_cast<pointer>(arr) ),
    bounds_( std::extent<ArrayType, I>::value... )
  { };  

  

  //============================================================================
  /// \name private data
  //============================================================================

  //! \brief a pointer to the data
  pointer data_ = nullptr;
  //! \brief the bounds of the array in multidimensional space
  bounds_type bounds_;
  
        
};
        
/// @}
    
}      // End namespace
}      // End namespace
