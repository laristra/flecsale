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
 * \file matrix.h
 * 
 * \brief Provides a dimensioned array which functions as a matrix.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/math/array/multi_array.h"
#include "ale/math/vector.h"

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
template <typename T, std::size_t D1, std::size_t D2> 
using row_major_matrix = array::multi_array< array::layouts::row_major, T, D1, D2>;


template <typename T, std::size_t D1, std::size_t D2> 
using col_major_matrix = array::multi_array< array::layouts::column_major, T, D1, D2>;


#if 0
////////////////////////////////////////////////////////////////////////////////
//! \brief Compute the dot product
//! \tparam T  The array base value type.
//! \tparam D  The array dimension.
//! \param[in] a  The first vector
//! \param[in] b  The other vector
//! \return The result of the operation
////////////////////////////////////////////////////////////////////////////////
template < 
  typename M, typename T, std::size_t D
>
std::enable_if_t< 
  std::is_same< typename M::layout_type, detail::row_major_layout >::value, M
> 
outer_product(const vector<T, D> &a, const vector<T, D> &b)
{

}

#endif

//template <typename T, std::size_t D>
//auto outer_product(const vector<T, D> &a, const vector<T, D> &b)
//{
//  row_major_matrix<T,D,D> tmp;
//  return tmp;
//}

//template <typename T, std::size_t D>
//auto outer_product(const vector<T, D> &a, const vector<T, D> &b)
//{
//  col_major_matrix<T,D,D> tmp;
//  return tmp;
//}

} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
