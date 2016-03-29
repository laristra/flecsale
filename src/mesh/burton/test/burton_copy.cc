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
 * \file burton_io.cc
 * 
 * \brief Tests io of the burton mesh.
 *
 ******************************************************************************/

//! user includes
#include "burton_test.h"

////////////////////////////////////////////////////////////////////////////////
//! \brief test copying the mesh
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, copy) {

  for(auto v : mesh_.vertices()){
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << std::endl;
  } // for

  auto mesh_copy = mesh_;

  for(auto v : mesh_copy.vertices()){
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << std::endl;
  } // for

  std::cout << CINCH_DUMP() << std::endl;

} // TEST

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
