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
 * \file qr.h
 * 
 * \brief Defines some functions for a qr solver
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "types.h"

namespace ale {
namespace linalg {

namespace detail {

///////////////////////////////////////////////////////////////////
/// Column Pivoting.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  typename T, 
  typename I,
  std::ptrdiff_t...MatArgs
>
auto get_next_col( const MatrixViewType<T, MatArgs...>&A, 
                   I row_pos,
                   I* p ) 
{
  
  // get the size type
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  auto max_loc = static_cast<size_type>(0);

  std::vector<T> col_norms(cols);
    
  // Compute the norms of the sub columns.
  for(size_type j = 0; j < cols; j++) {
    col_norms[j] = 0;

    for(size_type i = row_pos; i < rows; i++)
      col_norms[j] += A(i,p[j])*A(i,p[j]);
  }

  // Find the maximum location.
  auto max = static_cast<T>(0);
  for(size_type i = 0; i < cols; i++)
    if(col_norms[i] > max) {
      max = col_norms[i];
      max_loc = i;
    }

  // Collect garbge and return.
  return max_loc;
}

///////////////////////////////////////////////////////////////////
/// Determine the householder transformation.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  typename T, 
  typename I,
  std::ptrdiff_t...MatArgs
>
void householder( const MatrixViewType<T, MatArgs...> &A,
                  I row_pos, I col_pos, 
                  T * result)
{
  // get the size type
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  auto norm = static_cast<T>(0);

  for(size_type i = row_pos; i < rows; i++)
    norm += A(i,col_pos)*A(i,col_pos);

  if(norm == 0) return;

  norm = std::sqrt(norm);

  result[0] = A(row_pos,col_pos) - norm;

  for(size_type i = 1; i < (rows - row_pos); i++)
    result[i] = A(i+row_pos,col_pos);

  norm = 0;
  for(size_type i = 0; i < (rows - row_pos); i++)
    norm += result[i]*result[i];

  if(norm == 0) return;

  norm = std::sqrt(norm);

  for(size_type i = 0; i < (rows - row_pos); i++)
    result[i] *= (1.0/norm);
}


///////////////////////////////////////////////////////////////////
/// Apply the householder transformation.  VECTOR RHS VERION
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  template <typename, std::ptrdiff_t...> class VectorViewType,
  typename T, 
  typename I,
  std::ptrdiff_t...MatArgs,
  std::ptrdiff_t...VecArgs
>
void apply_householder( const MatrixViewType<T, MatArgs...> & A,
                        const VectorViewType<T, VecArgs...> & B, 
                        T * house, 
                        I row_pos, I *p )
{
    
  // get the size type
  using size_type = typename matrix_view<T>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  // Get the dimensions for the Q matrix.
  auto nn = rows - row_pos;

  // some temporary matrices
  std::vector<T> A_cpy( A.begin(), A.end() );
  std::vector<T> B_cpy( B.begin(), B.end() );
  std::vector<T> hhmat(nn*nn);

  // shape the matrix view
  auto A_cpy_view = MatrixViewType<T, MatArgs...>( A_cpy, A.bounds() );
      
  // Build the Q matrix from the Householder transform.
  for(size_type j = 0; j < nn; j++)
    for(size_type i = 0; i < nn; i++)
      if(i != j)
        hhmat[i+nn*j] = -2*house[j]*house[i];
      else
        hhmat[i+nn*j] = 1 - 2*house[j]*house[i];

  // Multiply by the Q matrix.
  for(size_type k = 0; k < cols; k++)
    for(size_type j = 0; j < nn; j++) {
        
      auto temp = static_cast<T>(0);

      for(size_type i = 0; i < nn; i++)
        temp += hhmat[i+nn*j]*A_cpy_view(i + row_pos,p[k]);
        
      A(j + row_pos,p[k]) = temp;
    }

  // Multiply the rhs by the Q matrix.
  for(size_type j = 0; j < nn; j++) {
        
    auto sum = static_cast<T>(0);
      
    for(size_type i = 0; i < nn; i++)
      sum += hhmat[i+nn*j]*B_cpy[i + row_pos];
        
    B[j + row_pos] = sum;
  }

}

///////////////////////////////////////////////////////////////////
/// Apply the householder transformation.  VECTOR RHS VERSION.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  template <typename, std::ptrdiff_t...> class VectorViewType,
  typename T, 
  typename I,
  std::ptrdiff_t...MatArgs,
  std::ptrdiff_t...VecArgs
>
void back_solve( const MatrixViewType<T, MatArgs...> & A,
                 const VectorViewType<T, VecArgs...> & B, 
                 I * p )
{

  // get epsilon
  constexpr auto eps = std::numeric_limits<T>::epsilon();

  // get the size type
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  size_type bottom = 0;

  // setup some temporary arrays
  std::vector<T> B_cpy( B.begin(), B.end() );
    
  // Find the first non-zero row from the bottom and start solving from here.
  for(size_type i = rows; i-- > 0;) {
    if( std::abs(A(i,p[cols - 1])) > eps ) {
      bottom = i;
      break;
    }
  }
    
  bottom = std::min( bottom, cols-1 );

  // Standard back solving routine starting at the first non-zero diagonal.
  for(size_type i = bottom+1; i-- > 0;) {
        
    auto sum = static_cast<T>(0);

    for(size_type j = cols; j-- > i+1;) 
      // if(j > i) 
      sum += B[p[j]]*A(i,p[j]);
      
    if ( std::abs(A(i,p[i])) > eps ) {
      auto temp = 1 / A(i,p[i]);
      B[ p[i] ] = (B_cpy[i] - sum) * temp;
    }
    else
      B[ p[i] ] = 0;
  }

}


} // namespace

///////////////////////////////////////////////////////////////////
/// computes the minimum-norm solution to a real linear least 
/// squares problem using a QR-based routine. VECTOR RHS VERSION.
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

  // get the size type
  using size_type = typename MatrixViewType<T,MatArgs...>::size_type;

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
  for(size_type i = 0; i < cols; i++) {

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


