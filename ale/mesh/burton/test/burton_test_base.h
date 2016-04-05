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
 * \file burton_test_base.h
 * 
 * \brief Defines a base test fixture.
 *
 ******************************************************************************/
#pragma once

//! system includes
#include <cinchtest.h>

//! user includes
#include "../../../mesh/burton/burton.h"

//! \brief the mesh type
using mesh_t   = ale::mesh::burton_mesh_t;

// some general using statements
using std::string;
using ale::mesh::write_mesh;
using ale::mesh::read_mesh;

////////////////////////////////////////////////////////////////////////////////
//! \brief base test fixture for burton
////////////////////////////////////////////////////////////////////////////////
class BurtonTestBase : public ::testing::Test {
protected:

  //---------------------------------------------------------------------------
  //! \brief get the output prefix for files
  //---------------------------------------------------------------------------
  auto output_prefix() 
  {
    auto test_info =
      ::testing::UnitTest::GetInstance()->current_test_info();
    string test_name = test_info->name();
    string case_name = test_info->test_case_name();
    return case_name + "." + test_name;
  }

};
