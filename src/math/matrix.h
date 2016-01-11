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
#include "ale/math/array.h"

namespace ale {
namespace math {

namespace detail {

////////////////////////////////////////////////////////////////////////////////
//!  \brief Define the layout of the matrix using a row-major convention
////////////////////////////////////////////////////////////////////////////////
struct row_major_layout {
  static constexpr 
  auto element( auto i, auto j, 
                auto /* size_i */, auto size_j ) noexcept
  { return i * size_j + j; }   
};


////////////////////////////////////////////////////////////////////////////////
//!  \brief Define the layout of the matrix using column major convention
////////////////////////////////////////////////////////////////////////////////
struct column_major_layout {
  static constexpr 
  auto element( auto i, auto j, 
                auto size_i, auto /* size_j */ ) noexcept
  { return j * size_i + i; }   
};

}


////////////////////////////////////////////////////////////////////////////////
//!  \brief The dimensioned_array type provides a general base for defining
//!  contiguous array types that have a specific dimension.
//!
//!  \tparam T The type of the array, e.g., P.O.D. type.
//!  \tparam D The dimension of the array, i.e., the number of elements
//!    to be stored in the array.
////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t D1, std::size_t D2> 
using row_major_matrix = math::array< detail::row_major_layout, T, D1, D2>;


template <typename T, std::size_t D1, std::size_t D2> 
using col_major_matrix = math::array< detail::column_major_layout, T, D1, D2>;


} // namespace
} // namespace

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
