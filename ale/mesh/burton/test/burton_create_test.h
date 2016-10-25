/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file burton_test.h
/// \brief Defines a test fixture.
////////////////////////////////////////////////////////////////////////////////

//! test include
#include "burton_test_base.h"


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh
////////////////////////////////////////////////////////////////////////////////
class burton_create : public burton_test_base {

public:


  ////////////////////////////////////////////////////////////////////////////////
  // type defines
  ////////////////////////////////////////////////////////////////////////////////

  //! \brief the mesh type
  using mesh_t = mesh_3d_t;
 
  //! \brief the mesh index type
  using size_t   = typename mesh_t::size_t;
  //! \brief the counter type
  using counter_t= typename mesh_t::counter_t;
 //! \brief the mesh float type
  using real_t   = typename mesh_t::real_t;
  //! \brief the mesh int type
  using integer_t= typename mesh_t::integer_t;
  //! \brief the mesh dimensions
  static constexpr auto num_dimensions = mesh_t::num_dimensions;

  //! \brief the vertex type
  using vertex_t = typename mesh_t::vertex_t;
  //! \brief the vertex type
  using face_t = typename mesh_t::face_t;

  ////////////////////////////////////////////////////////////////////////////////
  // the test parameters
  ////////////////////////////////////////////////////////////////////////////////

  //! \brief number of cells wide
  static constexpr size_t num_cells_x = 3;
  static constexpr size_t num_cells_y = 2;
  static constexpr size_t num_cells_z = 2;
  static constexpr auto num_cells = num_cells_x * num_cells_y * num_cells_z;


  //! \brief comput the number of vertices in each direction
  static constexpr auto num_vert_x = num_cells_x + 1;
  static constexpr auto num_vert_y = num_cells_y + 1;
  static constexpr auto num_vert_z = num_cells_z + 1;
  static constexpr auto num_vert = num_vert_x * num_vert_y * num_vert_z;

  // the number of faces in each direction
  static constexpr auto num_face_x = num_vert_x *num_cells_y*num_cells_z;
  static constexpr auto num_face_y = num_cells_x*num_vert_y *num_cells_z;
  static constexpr auto num_face_z = num_cells_x*num_cells_y*num_vert_z;
  //! \brief comput the number of faces in each direction
  static constexpr auto num_faces = num_face_x + num_face_y + num_face_z;

  //! \brief some test tolerance
  static constexpr real_t test_tolerance = ale::common::test_tolerance;


  // the vertex strides in each direction
  static constexpr auto vert_stride_x = 1;
  static constexpr auto vert_stride_y = vert_stride_x * num_vert_x;
  static constexpr auto vert_stride_z = vert_stride_y * num_vert_y;


  // a 3d to 1d indexing function for the vertices
  auto vert_index(size_t i, size_t j, size_t k) 
  {
    return i*vert_stride_x +  j*vert_stride_y +  k*vert_stride_z;
  };


  // the strides in each direction
  static constexpr auto face_x_stride_x = 1;
  static constexpr auto face_x_stride_y = face_x_stride_x * num_vert_x;
  static constexpr auto face_x_stride_z = face_x_stride_y * num_cells_y;

  static constexpr auto face_y_stride_x = 1;
  static constexpr auto face_y_stride_y = face_y_stride_x * num_cells_x;
  static constexpr auto face_y_stride_z = face_y_stride_y * num_vert_y;

  static constexpr auto face_z_stride_x = 1;
  static constexpr auto face_z_stride_y = face_z_stride_x * num_cells_x;
  static constexpr auto face_z_stride_z = face_z_stride_y * num_cells_y;

  // a 3d to 1d indexing function for the faces
  auto face_x_index(size_t i, size_t j, size_t k) 
  {
    return i*face_x_stride_x +  j*face_x_stride_y +  k*face_x_stride_z;
  };

  auto face_y_index(size_t i, size_t j, size_t k) 
  {
    return i*face_y_stride_x +  j*face_y_stride_y +  k*face_y_stride_z;
  };

  auto face_z_index(size_t i, size_t j, size_t k) 
  {
    return i*face_z_stride_x +  j*face_z_stride_y +  k*face_z_stride_z;
  };

};
