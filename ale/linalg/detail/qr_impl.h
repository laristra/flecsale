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

namespace ale {
namespace linalg {
namespace detail {

///////////////////////////////////////////////////////////////////
/// \brief Perform column pivoting.
///
/// \param [in,out] A  The system matrix.
/// \param [in]  row_pos  The starting row position.
/// \param [in]  p        The column pivots.
/// \return The maximum column pivot.
///
/// \tparam T  The value type.
/// \tparam I  The index type.
/// \tparam MatrixViewType  The type of the matrix view.
/// \tparam MatArgs  The deduced arguments of the matrix view.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  typename T, 
  typename I,
  typename J,
  std::ptrdiff_t...MatArgs
>
auto get_next_col( const MatrixViewType<T, MatArgs...>&A, 
                   I row_pos,
                   J* p ) 
{
  
  // get the size and coutner types
  using counter_type = typename MatrixViewType<T, MatArgs...>::counter_type;
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  auto max_loc = static_cast<size_type>(0);

  std::vector<T> col_norms(cols);
    
  // Compute the norms of the sub columns.
  for(counter_type j = 0; j < cols; j++) {
    col_norms[j] = 0;

    for(counter_type i = row_pos; i < rows; i++)
      col_norms[j] += A(i,p[j])*A(i,p[j]);
  }

  // Find the maximum location.
  auto max = static_cast<T>(0);
  for(counter_type i = 0; i < cols; i++)
    if(col_norms[i] > max) {
      max = col_norms[i];
      max_loc = i;
    }

  // Collect garbge and return.
  return max_loc;
}

///////////////////////////////////////////////////////////////////
/// \brief Determine the householder transformation to apply.
///
/// \param [in,out] A  The system matrix.
/// \param [in]  row_pos  The starting row position.
/// \param [in]  col_pos  The column pivots.
/// \param [out] result   The result of the transformation.
///
/// \tparam T  The value type.
/// \tparam I  The index type.
/// \tparam J  The index type.
/// \tparam MatrixViewType  The type of the matrix view.
/// \tparam MatArgs  The deduced arguments of the matrix view.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  typename T, 
  typename I,
  typename J,
  std::ptrdiff_t...MatArgs
>
void householder( const MatrixViewType<T, MatArgs...> &A,
                  I row_pos, J col_pos, 
                  T * result)
{
  // get the size and coutner types
  using counter_type = typename MatrixViewType<T, MatArgs...>::counter_type;
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  auto norm = static_cast<T>(0);

  for(counter_type i = row_pos; i < rows; i++)
    norm += A(i,col_pos)*A(i,col_pos);

  if(norm == 0) return;

  norm = std::sqrt(norm);

  result[0] = A(row_pos,col_pos) - norm;

  for(counter_type i = 1; i < (rows - row_pos); i++)
    result[i] = A(i+row_pos,col_pos);

  norm = 0;
  for(counter_type i = 0; i < (rows - row_pos); i++)
    norm += result[i]*result[i];

  if(norm == 0) return;

  norm = std::sqrt(norm);

  for(counter_type i = 0; i < (rows - row_pos); i++)
    result[i] *= (1.0/norm);
}


///////////////////////////////////////////////////////////////////
/// \brief Apply the householder transformation.
///
/// \param [in,out] A  The system matrix.
/// \param [in,out] B  The right hand side vector.
/// \param [in]  row_pos  The starting row position.
/// \param [in]  p        The column pivots.
/// \param [in]  house    The householder transformation to apply.
///
/// \tparam T  The value type.
/// \tparam I  The index type.
/// \tparam J  The index array type.
/// \tparam MatrixViewType  The type of the matrix view.
/// \tparam VectorViewType  The type of the vector view.
/// \tparam MatArgs,VecArgs  The deduced arguments of the matrix 
///                          and vector views.
///////////////////////////////////////////////////////////////////
template< 
  template <typename, std::ptrdiff_t...> class MatrixViewType,
  template <typename, std::ptrdiff_t...> class VectorViewType,
  typename T, 
  typename I,
  typename J,
  std::ptrdiff_t...MatArgs,
  std::ptrdiff_t...VecArgs
>
void apply_householder( const MatrixViewType<T, MatArgs...> & A,
                        const VectorViewType<T, VecArgs...> & B, 
                        T * house, 
                        I row_pos, J *p )
{
    
  // get the size and coutner types
  using counter_type = typename MatrixViewType<T, MatArgs...>::counter_type;
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

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
  for(counter_type j = 0; j < nn; j++)
    for(counter_type i = 0; i < nn; i++)
      if(i != j)
        hhmat[i+nn*j] = -2*house[j]*house[i];
      else
        hhmat[i+nn*j] = 1 - 2*house[j]*house[i];

  // Multiply by the Q matrix.
  for(counter_type k = 0; k < cols; k++)
    for(counter_type j = 0; j < nn; j++) {
        
      auto temp = static_cast<T>(0);

      for(counter_type i = 0; i < nn; i++)
        temp += hhmat[i+nn*j]*A_cpy_view(i + row_pos,p[k]);
        
      A(j + row_pos,p[k]) = temp;
    }

  // Multiply the rhs by the Q matrix.
  for(counter_type j = 0; j < nn; j++) {
        
    auto sum = static_cast<T>(0);
      
    for(counter_type i = 0; i < nn; i++)
      sum += hhmat[i+nn*j]*B_cpy[i + row_pos];
        
    B[j + row_pos] = sum;
  }

}

///////////////////////////////////////////////////////////////////
/// \brief Apply back substitution to get the solution.
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

  // get the size and coutner types
  using counter_type = typename MatrixViewType<T, MatArgs...>::counter_type;
  using size_type = typename MatrixViewType<T, MatArgs...>::size_type;

  // some matrix dimensions
  auto rows = A.template extent<0>();
  auto cols = A.template extent<1>();

  size_type bottom = 0;

  // setup some temporary arrays
  std::vector<T> B_cpy( B.begin(), B.end() );
    
  // Find the first non-zero row from the bottom and start solving from here.
  for(counter_type i = rows; i-- > 0;) {
    if( std::abs(A(i,p[cols - 1])) > eps ) {
      bottom = i;
      break;
    }
  }
    
  bottom = std::min( bottom, cols-1 );

  // Standard back solving routine starting at the first non-zero diagonal.
  for(counter_type i = bottom+1; i-- > 0;) {
        
    auto sum = static_cast<T>(0);

    for(counter_type j = cols; j-- > i+1;) 
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
} // namespace
} // namespace


