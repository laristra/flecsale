/////////////////////////////////////////////////////////////////////
/// \author Marc R.J. "T-Bone"  Charest
///
/// Functions for solving linear systems based on QR factorization.
///
/// \date Thursday, August 12 2010
/////////////////////////////////////////////////////////////////////
#pragma once

namespace ale {
namespace math {

namespace detail {

///////////////////////////////////////////////////////////////////
/// Column Pivoting.
///////////////////////////////////////////////////////////////////
template< typename T, typename I >
auto get_next_col( const utils::multi_array_ref<T> &A, 
                   I row_pos,
                   I* p ) 
{
  
  // get the size type
  using size_type = utils::multi_array_ref<T>::size_type;

  // some matrix dimensions
  auto rows = A.dimension(0);
  auto cols = A.dimension(1);

  auto max_loc = static_cast<size_type>(0);

  vector<T> col_norms(cols);
    
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
template< typename T, typename I >
void householder( const utils::multi_array_ref<T> &A,
                  I row_pos, I col_pos, 
                  T * result)
{
  // get the size type
  using size_type = utils::multi_array_ref<T>::size_type;

  // some matrix dimensions
  auto rows = A.dimension(0);
  auto cols = A.dimension(1);

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

   norm = sqrt(norm);

  for(size_type i = 0; i < (rows - row_pos); i++)
    result[i] *= (1.0/norm);
}


///////////////////////////////////////////////////////////////////
/// Apply the householder transformation.  VECTOR RHS VERION
///////////////////////////////////////////////////////////////////
template< typename T, typename I >
void apply_householder( const utils::multi_array_ref<T> &A,
                        T * B, 
                        T * house, 
                        I row_pos, I *p )
{
    
  // get the size type
  using size_type = utils::multi_array_ref<T>::size_type;

  // some matrix dimensions
  auto rows = A.dimension(0);
  auto cols = A.dimension(1);

  // Get the dimensions for the Q matrix.
  auto nn = rows - row_pos;

  // some temporary matrices
  math::Matrix2D A_cpy(A);
  std::vector<T> B_cpy( B, B+rows );
  std::vector<T> hhmat(nn*nn);
      
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
        temp += hhmat[i+nn*j]*A_cpy(i + row_pos,p[k]);
        
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
template< typename T, typename I >
void back_solve( const utils::multi_array_ref<T> &A,
                 T * B, 
                 I * p )
{

  // get epsilon
  constexpr auto eps = std::numeric_limits<T>::epsilon();

  // get the size type
  using size_type = utils::multi_array_ref<T>::size_type;

  // some matrix dimensions
  auto rows = A.dimension(0);
  auto cols = A.dimension(1);

  size_type bottom;

  // setup some temporary arrays
  std::vector<T> B_cpy( B, B+rows );
    
  // Find the first non-zero row from the bottom and start solving from here.
  for(size_type i = rows - 1; i >= 0; i--) {
    if( std::abs(A(i,p[cols - 1])) > eps ) {
      bottom = i;
      break;
    }
  }
    
  bottom = std::min( bottom, cols-1 );

  // Standard back solving routine starting at the first non-zero diagonal.
  for(size_type i = bottom; i >= 0; i--) {
        
    auto sum = static_cast<T>(0);

    for(size_type j = cols - 1; j >= 0; j--)
      if(j > i)
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
template< typename T >
solve_qr( utils::mutli_array_ref<T> &A, 
          utils::array_ref<T> &B) 
{

  auto rows = A.dimension(0);
  auto cols = A.dimension(1);
    
  // initial checks
  if (rows < 1 || cols < 1) 
    raise_runtime_error("Incorect matrix sizes");

  if ( B.size() != rows ) 
    raise_runtime_error("RHS vector wrong size");

  // householder vector
  std::vector<T> v(rows);


  // Initial permutation vector.
  auto jpvt = vector<size_type>(cols);
  std::iota( jpvt.begin(), jpvt.end(), static_cast<size_type>(0) );
  
  // Apply rotators to make R and Q'*b 
  for(size_type i = 0; i < cols; i++) {

    // work[0:cols)   -> col_norms
    auto max_loc = get_next_col(A, i, jpvt.data());

    if(max_loc >= 0)
      std::swap(jpvt[i], jpvt[max_loc]);

    // work -> Q matrix, A copy, B copy
    detail::householder(A, i, jpvt[i], v.data());
    detail::apply_householder(A, B.data(), v.data(), i, jpvt.data());

  }


  // Back solve Rx = Q'*b 
  detail::back_solve(A, B, jpvt);
}


} // namespace
} // namespace


