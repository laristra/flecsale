/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Defines some functions for a qr solver.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "detail/qr_impl.h"

#include "types.h"

namespace ale {
namespace linalg {


///////////////////////////////////////////////////////////////////
/// \brief Computes the minimum-norm solution to a real linear least 
/// squares problem using a QR-based routine.
///
/// Solves for `x` in `A x = B`
///
/// \param [in,out] A  The system matrix.
/// \param [in,out] B  On entry, the right hand side vector.  On 
///                    exit, the solution vector.
///
/// \tparam T  The value type.
/// \tparam MatrixViewType  The type of the matrix view.
/// \tparam VectorViewType  The type of the vector view.
/// \tparam MatArgs,VecArgs  The deduced arguments of the matrix 
///                          and vector views.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  template <typename, std::ptrdiff_t...> class VectorViewType,
  typename T, std::ptrdiff_t...MatArgs, std::ptrdiff_t...VecArgs
>
void qr ( 
  MatrixViewType<T,MatArgs...> A, 
  VectorViewType<T,VecArgs...> B 
) {

  // some static assertions
  static_assert( A.rank() == 2, "System matrix must have rank 2" );
  static_assert( B.rank() == 1, "Right-hand-side vector must have rank 1" );

  // get the size and coutner types
  using counter_type = typename MatrixViewType<T, MatArgs...>::counter_type;
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // the dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();
    
  // initial checks
  if (rows < 1 || cols < 1) 
    raise_runtime_error("Incorect matrix sizes");

  if ( B.template extent<0>() != rows ) 
    raise_runtime_error("RHS vector wrong size");

  // householder vector
  std::vector<T> v(rows);


  // Initial permutation vector.
  auto jpvt = std::vector<size_type>(cols);
  std::iota( jpvt.begin(), jpvt.end(), static_cast<size_type>(0) );
  
  // Apply rotators to make R and Q'*b 
  for(counter_type i = 0; i < cols; i++) {

    // work[0:cols)   -> col_norms
    auto max_loc = detail::get_next_col(A, i, jpvt.data());

    if(max_loc >= 0)
      std::swap(jpvt[i], jpvt[max_loc]);

    // work -> Q matrix, A copy, B copy
    detail::householder(A, i, jpvt[i], v.data());
    detail::apply_householder(A, B, v.data(), i, jpvt.data());

  }


  // Back solve Rx = Q'*b 
  detail::back_solve(A, B, jpvt.data());
}


} // namespace
} // namespace


