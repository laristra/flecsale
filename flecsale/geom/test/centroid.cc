/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Tests related to the centroid operation.
///
////////////////////////////////////////////////////////////////////////////////

// preprocessor directive for exception handling
#define ENABLE_EXCEPTIONS

// system includes
#include<cinchtest.h>
#include<vector>

// user includes
#include "flecsale/common/types.h"
#include "flecsale/geom/centroid.h"
#include "flecsale/geom/point.h"

// explicitly use some stuff
using std::cout;
using std::endl;
using std::vector;

using namespace flecsale;
using namespace flecsale::geom;

//! the real type
using real_t = common::real_t;
// ! the test tolerance
using common::test_tolerance;

//! the point type
using point_2d_t = point<real_t, 2>;

//! the error types used
using utils::ExceptionNotImplemented;

//=============================================================================
//! \brief Test centroid operator.
//=============================================================================
TEST(centroid, 2d) 
{

  vector<point_2d_t> points = { {0, 0}, {1, 0}, {1, 1}, {0, 1} };
  
  auto cx1 = centroid( 
    point_2d_t{0, 0}, point_2d_t{1, 0}, point_2d_t{1, 1}, point_2d_t{0, 1} 
  );
  auto cx2 = centroid( points );

  ASSERT_NEAR( 0.5, cx1[0], test_tolerance ) << " Centroid calculation wrong ";
  ASSERT_NEAR( 0.5, cx2[0], test_tolerance ) << " Centroid calculation wrong ";

  ASSERT_NEAR( 0.5, cx1[1], test_tolerance ) << " Centroid calculation wrong ";
  ASSERT_NEAR( 0.5, cx2[1], test_tolerance ) << " Centroid calculation wrong ";

} // TEST

/* Test that centroid operator correctly handles unsupported
 * input vector dimensionalities */
using point_3d_t = point<real_t, 3>;
using point_5d_t = point<real_t, 5>;

TEST(centroid, HandlesUnsupportedDims)
{
  vector<point_3d_t> points_3d = { {0, 0, 0}, {1, 0, 7} };
  vector<point_5d_t> points_5d = { {2, 0, 9, 0, 3}, {8 ,7, 9, 0, 7} };

  ASSERT_THROW( centroid( points_3d ), ExceptionNotImplemented);
  ASSERT_THROW( centroid( points_5d ), ExceptionNotImplemented);
} // TEST
