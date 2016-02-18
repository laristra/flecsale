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
#include "burton_test.h"

// using statements
using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////////////
//! \brief dump the mesh to std out
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, dump) {
  mesh_.dump();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh connectivity
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, mesh) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Vertices in mesh:" << endl;
  for(auto v : mesh_.vertices()){
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;
  }

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Edges in mesh:" << endl;
  for(auto e : mesh_.edges()){
    CINCH_CAPTURE() << "----------- edge id: " << e.id()
                    << " with midpoint " << e->midpoint() << endl;
  }

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Corners in mesh:" << endl;

  for(auto c : mesh_.corners()) {
    CINCH_CAPTURE() << "----------- corner id: " << c.id() << endl;

    for(auto e: mesh_.edges(c)) {
      CINCH_CAPTURE() << "      ++++ edge id: " << e.id() << endl;
    } // for
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each cell:" << endl;
  for(auto c : mesh_.cells()) {
    CINCH_CAPTURE() << "-----------Edges for cell id: " << c.id()
                    << " with centroid " << c->centroid() << endl;
  }

  //CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED("burton.blessed"));

  cout << CINCH_DUMP() << endl;

} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test the accessors
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, accessors) {
  register_state(mesh_, "pressure", cells, real_t, persistent);
  register_state(mesh_, "density", cells, real_t);
  register_state(mesh_, "total energy", cells, real_t, persistent);
  register_state(mesh_, "velocity", vertices, vector_t, persistent);
  register_state(mesh_, "H", edges, vector_t);

  cout << "Accessing state with type real_t:" << endl;

  auto vr = access_type(mesh_, real_t);
  for(auto v: vr) {
    cout << "\t" << v.label() << " has type real_t" << endl;
  } // for

  cout << endl;

  cout << "Accessing state with type real_t at cells:" << endl;

  auto va = access_type_if(mesh_, real_t, is_at(cells));
  for(auto v: va) {
    cout << "\t" << v.label() <<
      " has type real_t and is at cells" << endl;
  } // for

  cout << endl;

  cout << "Accessing persistent state with type real_t at cells:" <<
    endl;

  auto vp = access_type_if(mesh_, real_t, is_persistent_at(cells));

  for(auto v: vp) {
    cout << "\t" << v.label() <<
      " has type real_t and is persistent at cells" << endl;
  } // for

  cout << endl;

  cout << "Accessing state with type vector_t at vertices:" << endl;

  auto vv = access_type_if(mesh_, vector_t, is_at(vertices));
  for(auto v: vv) {
    cout << "\t" << v.label() <<
      " has type vector_t and is at vertices" << endl;
  } // for

  cout << endl;

  cout << "Accessing persistent state with type vector_t at vertices:"
            << endl;

  auto vpv = access_type_if(mesh_, vector_t, is_persistent_at(vertices));
  for(auto v: vpv) {
    cout << "\t" << v.label() <<
      " has type vector_t and is persistent at vertices" << endl;
  } // for

  cout << endl;

  cout << "Accessing state with type vector_t at edges:" << endl;

  auto ve = access_type_if(mesh_, vector_t, is_at(edges));
  for(auto v: ve) {
    cout << "\t" << v.label() <<
      " has type vector_t and is at edges" << endl;
  } // for

  cout << endl;

  cout << "Accessing state with type point_t at vertices:" << endl;

  auto pv = access_type_if(mesh_, point_t, is_at(vertices));
  for(auto v: pv) {
    cout << "\t" << v.label() <<
      " has type point_t and is at vertices" << endl;
  } // for

  cout << endl;

  cout << "Accessing persistent state with type point_t at vertices:"
            << endl;

  auto ppv = access_type_if(mesh_, point_t, is_persistent_at(vertices));
  for(auto v: ppv) {
    cout << "\t" << v.label() <<
      " has type point_t and is persistent at vertices" << endl;
  } // for

  cout << endl;

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the state
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, state) {

  register_state(mesh_, "pressure", cells, real_t, persistent);
  register_state(mesh_, "velocity", vertices, vector_t, persistent);
  register_state(mesh_, "H", edges, vector_t);
//  register_state(mesh_, "cornerdata", corners, int32_t);
//  register_state(mesh_, "wedgedata", wedges, bool);

  auto p = access_state(mesh_, "pressure", real_t);
  auto velocity = access_state(mesh_, "velocity", vector_t);
  auto H = access_state(mesh_, "H", vector_t);
//  auto cd = access_state(mesh_, "cornerdata", int32_t);
//  auto wd = access_state(mesh_, "wedgedata", bool);

  // cells
  ASSERT_EQ(4, mesh_.num_cells());
  for(auto c: mesh_.cells()) {
    p[c] = c.id();
  } // for

  for(auto c: mesh_.cells()) {
    ASSERT_EQ(c.id(), p[c]);
  } // for

  // vertices
  ASSERT_EQ(9, mesh_.num_vertices());
  for (auto v: mesh_.vertices()) {
    velocity[v][0] = v.id();
    velocity[v][1] = 2.0*v.id();
  } // for

  for (auto v: mesh_.vertices()) {
    ASSERT_EQ(v.id(), velocity[v][0]);
    ASSERT_EQ(2.0*v.id(), velocity[v][1]);
  } // for

  // edges
  ASSERT_EQ(12, mesh_.num_edges());
  for (auto e: mesh_.edges()) {
    H[e][0] = e.id()*e.id();
    H[e][1] = e.id()*e.id()*e.id();
  } // for

  for (auto e: mesh_.edges()) {
    ASSERT_EQ(e.id()*e.id(), H[e][0]);
    ASSERT_EQ(e.id()*e.id()*e.id(), H[e][1]);
  } // for

} // TEST_F





/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
