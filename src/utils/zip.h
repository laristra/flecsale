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
 * \file zip.h
 * 
 * \brief Provide a zip-like iterator for range-based fors.
 *
 ******************************************************************************/
#pragma once

//#define USE_BOOST_ZIP

// system includes
#ifdef USE_BOOST_ZIP
#  include <boost/iterator/zip_iterator.hpp>
#  include <boost/range.hpp>
#else
#  include <iterator>
#  include <tuple>
#  include <utility>
#endif

namespace ale {
namespace utils {

#ifdef USE_BOOST_ZIP

////////////////////////////////////////////////////////////////////////////////
//! \brief Combine iterators together using boost
////////////////////////////////////////////////////////////////////////////////

using boost::get;

template <typename... T>
auto zip(T&&... containers)
{
  // get the number of arguments
  auto n = sizeof...(T);

  // check all the sizes
  size_t sizes[] = { containers.size()... };
  auto first_size = sizes[1];
  for ( size_t i=2; i<n; i++ )
    assert( first_size == sizes[i] && " size mismatch" ); 

  // do the zip
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}


#else 

////////////////////////////////////////////////////////////////////////////////
//! \brief Combine iterators together
//!
//! Got this from 
//! http://codereview.stackexchange.com/questions/30846/zip-like-functionality-with-c11s-range-based-for-loop
//!
////////////////////////////////////////////////////////////////////////////////

// need this in the open
using std::get;

// helper for tuple_subset and tuple_tail (from http://stackoverflow.com/questions/8569567/get-part-of-stdtuple)
template <size_t... n>
struct ct_integers_list {
  template <size_t m>
  struct push_back
  {
    using type = ct_integers_list<n..., m>;
  };
};

template <size_t max>
struct ct_iota_1
{
  using type = typename ct_iota_1<max-1>::type::template push_back<max>::type;
};

template <>
struct ct_iota_1<0>
{
  using type = ct_integers_list<>;
};

// return a subset of a tuple
template <size_t... indices, typename Tuple>
auto tuple_subset(Tuple&& tpl, ct_integers_list<indices...>)
{
  return std::make_tuple(std::get<indices>(std::forward<Tuple>(tpl))...);
  // this means:
  //   make_tuple(get<indices[0]>(tpl), get<indices[1]>(tpl), ...)
}

// return the tail of a tuple
template <typename Head, typename... Tail>
auto tuple_tail(const std::tuple<Head, Tail...>& tpl)
{
  return tuple_subset(tpl, typename ct_iota_1<sizeof...(Tail)>::type());
  // this means:
  //   tuple_subset<1, 2, 3, ..., sizeof...(Tail)-1>(tpl, ..)
}

// increment every element in a tuple (that is referenced)
template<std::size_t I = 0, typename... Tp>
typename std::enable_if<I == sizeof...(Tp), void>::type
increment(std::tuple<Tp...>& t)
{ }

template<std::size_t I = 0, typename... Tp>
typename std::enable_if<(I < sizeof...(Tp)), void>::type
increment(std::tuple<Tp...>& t)
{
  std::get<I>(t)++ ;
  increment<I + 1, Tp...>(t);
}

// check equality of a tuple
template<typename T1>
bool not_equal_tuples( const std::tuple<T1>& t1,  const std::tuple<T1>& t2 )
{
  return (std::get<0>(t1) != std::get<0>(t2));
}

template<typename T1, typename... Ts>
bool not_equal_tuples( const std::tuple<T1, Ts...>& t1,  const std::tuple<T1, Ts...>& t2 )
{
  return (std::get<0>(t1) != std::get<0>(t2)) && not_equal_tuples( tuple_tail(t1), tuple_tail(t2) );
}

// check equality of a tuple
template<typename T1>
bool equal_tuples( const std::tuple<T1>& t1,  const std::tuple<T1>& t2 )
{
  return (std::get<0>(t1) == std::get<0>(t2));
}

template<typename T1, typename... Ts>
bool equal_tuples( const std::tuple<T1, Ts...>& t1,  const std::tuple<T1, Ts...>& t2 )
{
  return (std::get<0>(t1) == std::get<0>(t2)) && equal_tuples( tuple_tail(t1), tuple_tail(t2) );
}

// dereference a subset of elements of a tuple (dereferencing the iterators)
template <size_t... indices, typename Tuple>
auto dereference_subset(Tuple& tpl, ct_integers_list<indices...>)
{
  return std::tie(*std::get<indices-1>(tpl)...);
}

// dereference every element of a tuple (applying operator* to each element, and returning the tuple)
template<typename... Ts>
auto dereference_tuple(std::tuple<Ts...>& t1)
{
  return dereference_subset( t1, typename ct_iota_1<sizeof...(Ts)>::type());
}


// a reference remover helper template
// from cppreference.com: 
template< class T >
using remove_reference_t = typename std::remove_reference<T>::type;


// a struct to get the value type of a referenced or unreferenced object
template <class T, class Enable = void>
struct value_type
{ };

template <class T>
struct value_type< T, 
                   typename std::enable_if< !std::is_reference<T>::value >::type >
{
  // use this to use by value, but i dont like that
  // using type = typename T::value_type;

  // always use references
  using type = typename T::value_type&;
};


template <class T>
struct value_type< T, 
                   typename std::enable_if< std::is_reference<T>::value >::type >
{
  using type =  typename remove_reference_t<T>::value_type&;
};

template< class T >
using value_type_t = typename value_type<T>::type;


// the zipper class with iterator
template< typename T1, typename... Ts >
class zipper
{
public:

  class iterator : std::iterator< std::forward_iterator_tag, 
                                  std::tuple<value_type_t<T1>, value_type_t<Ts>...> >
  {
  protected:
    std::tuple<typename remove_reference_t<T1>::iterator, 
               typename remove_reference_t<Ts>::iterator...> current;
  public:

    explicit iterator(  typename remove_reference_t<T1>::iterator s1, 
                        typename remove_reference_t<Ts>::iterator... s2 ) : 
      current(s1, s2...) {};

    iterator( const iterator& rhs ) :  current(rhs.current) {};

    auto & operator++() {
      increment(current);
      return *this;
    }

    auto operator++(int) {
      auto a = *this;
      increment(current);
      return a;
    }

    bool operator!=( const iterator& rhs ) {
      return not_equal_tuples(current, rhs.current);
    }

    bool operator==( const iterator& rhs ) {
      return equal_tuples(current, rhs.current);
    }

    typename iterator::value_type operator*() {
      return dereference_tuple(current);
    }
  };


  zipper( T1& a, Ts&... b):
    begin_( a.begin(), (b.begin())...), 
    end_( a.end(), (b.end())...) {};

  zipper(const zipper<T1, Ts...>& a) :
    begin_( a.begin_ ), 
    end_( a.end_ ) {};

  template<typename U1, typename... Us>
  auto & operator=(const zipper<U1, Us...>& rhs) {
    begin_ = rhs.begin_;
    end_ = rhs.end_;
    return *this;
  }

  auto& begin() 
  { return begin_; }

  auto& end() 
  { return end_; }



  iterator begin_;
  iterator end_;
};


//from cppreference.com: 
template <class T>
struct special_decay
{
  using type = typename std::decay<T>::type;
};

//allows the use of references:

template <class T>
struct special_decay<std::reference_wrapper<T>>
{
  using type = T&;
};

template <class T>
using special_decay_t = typename special_decay<T>::type;

//allows template type deduction for zipper:
template <class... Types>
auto zip(Types&&... args)
{

  // get the number of arguments
  auto n = sizeof...(Types);

  // check all the sizes
  size_t sizes[] = { args.size()... };
  auto first_size = sizes[1];
  for ( size_t i=2; i<n; i++ )
    assert( first_size == sizes[i] && " size mismatch" ); 

  return zipper<special_decay_t<Types>...>(std::forward<Types>(args)...);
}

#endif

} // namespace
} // namespace
