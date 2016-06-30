/*~-------------------------------------------------------------------------~~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/    /_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
/*!
 *
 * \file fixed_vector.cc
 * 
 * \brief Tests related to a fixed vector.
 *
 ******************************************************************************/

// system includes
#include<cinchtest.h>
#include<iostream>

// user includes
#include "ale/utils/array_view.h"

using namespace ale::utils;

///////////////////////////////////////////////////////////////////////////////
//! \brief Test the construction of bounds_t's.
///////////////////////////////////////////////////////////////////////////////
TEST(array_view, bounds) 
{

  bounds_t< 1 > bounds_1d_var(1);
  bounds_t< 1 > bounds_1d_list( {1} );

  bounds_t< 2 > bounds_2d_var(1, 2);
  bounds_t< 2 > bounds_2d_list( {1, 2} );
  bounds_t< 2 > bounds_2d_copy( bounds_2d_list );

  ASSERT_EQ( bounds_2d_var, bounds_2d_list );

  bounds_2d_copy[0] = 2;
  ASSERT_NE( bounds_2d_list, bounds_2d_copy );

} // TEST


///////////////////////////////////////////////////////////////////////////////
//! \brief Test the construction of index_t's.
///////////////////////////////////////////////////////////////////////////////
TEST(array_view, index) 
{
  index_t< 1 > index_1d_var(1);
  index_t< 1 > index_1d_list( {1} );
  index_t< 1 > index_1d_copy( index_1d_list );

  index_t< 2 > index_2d_var(1, 2);
  index_t< 2 > index_2d_list( {1, 2} );
  index_t< 2 > index_2d_copy( index_2d_list );
  
  index_t< 2 > index_2d_agg1{1, 2};
  index_t< 2 > index_2d_agg2 = {1, 2};
  
  ASSERT_EQ( index_2d_var, index_2d_list );

  index_2d_copy[0] = 2;
  ASSERT_NE( index_2d_list, index_2d_copy );

  index_1d_var++;
  ASSERT_EQ( index_1d_var[0], 2 );

} // TEST

///////////////////////////////////////////////////////////////////////////////
//! \brief Test the construction of index_t's.
///////////////////////////////////////////////////////////////////////////////
TEST(array_view, arithmetic) 
{

  auto bnd1 = index_t<3>{ 3, 1, 4 };
  auto idx = index_t<3>{ 2, -1, 0 };
  auto  bnd2 = bnd1 + idx;  // bnd2 is { 5, 0, 4 }
  ASSERT_EQ( bnd2, index_t<3>({ 5, 0, 4 }) );
  bnd1 -= idx;  // bnd1 is { 1, 2, 4 }
  ASSERT_EQ( bnd1, index_t<3>({ 1, 2, 4 }) );

  index_t<2> idx2{ 2, 3 };
  auto res = idx2 * 1.5;  // res is {3, 4}
  ASSERT_EQ( res, index_t<2>({ 3, 4 }) );
 
} // TEST

///////////////////////////////////////////////////////////////////////////////
//! \brief Test the construction of bounds_iterators's.
///////////////////////////////////////////////////////////////////////////////
TEST(array_view, bounds_iterator) 
{

  auto bnds = bounds_t<3>{ 3, 1, 4 };
  auto it = bnds.begin();
  ASSERT_TRUE( bnds.begin() != bnds.end() );
} // TEST



///////////////////////////////////////////////////////////////////////////////
//! \brief Test the construction of array_views's.
///////////////////////////////////////////////////////////////////////////////
TEST(array_view, constructor) 
{
  auto vec = std::vector<int>(10);
  auto view = array_view<int>{ vec };
  view[0] = 42;
  int v = vec[0]; // v == 42
  ASSERT_EQ( vec[0], 42 );

  // create a multidimensional array
  int r[3][1][2];
  for ( int i=0; i<3; i++ ) 
    for ( int j=0; j<1; j++ )
      for ( int k=0; k<2; k++ ) 
        r[i][j][k] = i + j + k;
        
  array_view<int, 3> mav{ r };  // av.bounds() is {3, 1, 2}
  array_view<int> fav{ mav };  // flatten multi-dimensional array
  
  // check ordering
  int cnt = 0;
  for ( int i=0; i<3; i++ )
    for ( int j=0; j<1; j++ )
      for ( int k=0; k<2; k++ ) {
        auto a = r[i][j][k];
        auto b = mav[{i,j,k}];
        auto c = fav[cnt++];
        auto d = mav(i, j, k);
        ASSERT_EQ( a, b );
        ASSERT_EQ( a, c );
        ASSERT_EQ( a, d );
      }
  
  // check copy constructor
  array_view<const int, 3> c_mav(mav);

  auto view_w_bounds_1 = array_view<int,2>( {2,5}, vec );
  auto view_w_bounds_2 = array_view<int,2>( {2,5}, vec.data() );

  array_view<const int,2> const_view;
  const_view = view_w_bounds_1;

} // TEST


///////////////////////////////////////////////////////////////////////////////
//! \brief Test the slicing array_views's.
///////////////////////////////////////////////////////////////////////////////
TEST(array_view, slice) 
{

  // create a multidimensional array
  int r[3][1][2];
  for ( int i=0; i<3; i++ ) 
    for ( int j=0; j<1; j++ )
      for ( int k=0; k<2; k++ ) 
        r[i][j][k] = i + j + k;
        
  array_view<int, 3> mav{ r };  // av.bounds() is {3, 1, 2}

  auto slice = mav.at( 0 );

  // check ordering
  for ( int j=0; j<1; j++ )
    for ( int k=0; k<2; k++ ) {
      auto a = r[0][j][k];
      auto b = slice[{j,k}];
      ASSERT_EQ( a, b );
    }

  auto orig = index_t<3>({1,0,1});
  auto sec = mav.section( orig );

  for ( int i=0; i<2; i++ ) 
    for ( int j=0; j<1; j++ )
      for ( int k=0; k<2; k++ ) {
        auto new_idx = orig + index_t<3>{i, j, k};
        auto a = mav[ new_idx ];
        auto b = sec[{i,j,k}];
        std::cout << "sec1 " << a << " " << b << std::endl;
        ASSERT_EQ( a, b );
      }

  auto sec2 = mav.section( orig, {2,1,2} );

  for ( int i=0; i<2; i++ ) 
    for ( int j=0; j<1; j++ )
      for ( int k=0; k<2; k++ ) {
        auto new_idx = orig + index_t<3>{i, j, k};
        auto a = mav[ new_idx ];
        auto b = sec2[{i,j,k}];
        std::cout << "sec2 " << a << " " << b << std::endl;
        ASSERT_EQ( a, b );
      }

} // TEST


/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
