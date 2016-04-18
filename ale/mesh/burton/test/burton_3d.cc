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
 * \file burton.cc
 * 
 * \brief Tests general features of the burton mesh.
 *
 ******************************************************************************/


//! user includes
#include "burton_nd_test.h"

// using statements
using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////////////
//! \brief dump the mesh to std out
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_3d, dump) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "dump" << endl;
  CINCH_CAPTURE() << separator;

  mesh_.dump();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh connectivity
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_3d, mesh) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "mesh" << endl;

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Vertices in mesh:" << endl;

  for(auto v : mesh_.vertices())
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Edges in mesh:" << endl;

  for(auto e : mesh_.edges())
    CINCH_CAPTURE() << "----------- edge id: " << e.id()
                    << " with midpoint " << e->midpoint() << endl;

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Faces in mesh:" << endl;

  for(auto f : mesh_.faces()) 
    CINCH_CAPTURE() << "----------- faces id: " << f.id()
                    << " with centroid " << f->centroid() << endl;

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Cells in mesh:" << endl;

  for(auto c : mesh_.cells()) 
    CINCH_CAPTURE() << "----------- cell id: " << c.id()
                    << " with centroid " << c->centroid() << endl;

  cout << CINCH_DUMP() << endl;

} // TEST_F



////////////////////////////////////////////////////////////////////////////////
//! \brief A final test to compare the blessed file and do CINCH_DUMP().
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_3d, cinch_dump) {
  cout << CINCH_DUMP() << endl;
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED("burton_3d.blessed"));
}

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
