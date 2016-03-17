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
TEST_F(Burton, mesh) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "mesh" << endl;
  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Vertices in mesh:" << endl;
  for(auto v : mesh_.vertices()){
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Edges in mesh:" << endl;
  for(auto e : mesh_.edges()){
    CINCH_CAPTURE() << "----------- edge id: " << e.id()
      << " with midpoint " << e->midpoint() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Corners in mesh:" << endl;

  for(auto c : mesh_.corners()) {
    CINCH_CAPTURE() << "----------- corner id: " << c.id() << endl;

    for(auto e: mesh_.edges(c)) {
      CINCH_CAPTURE() << "       ++++ edge id: " << e.id() << endl;
    } // for
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Wedges in mesh:" << endl;

  for(auto w : mesh_.wedges()) {
    CINCH_CAPTURE() << "----------- wedge id: " << w.id() << endl;

    for(auto c: mesh_.cells(w)) {
      CINCH_CAPTURE() << "       ++++ cell id: " << c.id() << endl;
    } // for

    for(auto e: mesh_.edges(w)) {
      CINCH_CAPTURE() << "       ++++ edge id: " << e.id() << endl;
    } // for

    for(auto v: mesh_.vertices(w)) {
      CINCH_CAPTURE() << "       ++++ vertex id: " << v.id() << endl;
    } // for
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each vertex:" << endl;
  for(auto v: mesh_.vertices()) {
    CINCH_CAPTURE() << "^^^^^^^^Vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;
    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w: mesh_.wedges(v)) {
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;
    } // for
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each cell:" << endl;
  for(auto c : mesh_.cells()) {
    CINCH_CAPTURE() << "^^^^^^^^Cell id: " << c.id()
      << " with centroid " << c->centroid() << endl;

    CINCH_CAPTURE() << "    ----Edges:" << endl;
    for(auto e : mesh_.edges(c)){
      CINCH_CAPTURE() << "    ++++ edge id: " << e.id()
        << " with midpoint " << e->midpoint() << endl;
    } // for

    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w : mesh_.wedges(c)){
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;
    } // for

    CINCH_CAPTURE() << "    ----Corners:" << endl;
    for(auto cnr : mesh_.corners(c)){
      CINCH_CAPTURE() << "    ++++ corner id: " << cnr.id() << endl;
    } // for
  } // for

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh geometry functions
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, geometry) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "geometry" << endl;
  CINCH_CAPTURE() << separator;
  for(auto c: mesh_.cells()) {
    auto xc = c->centroid();
    CINCH_CAPTURE() << "---- cell id: " << c.id()
      << " with centroid " << xc << endl;
    for(auto v : mesh_.vertices(c)){
      auto xv = v->coordinates();
      CINCH_CAPTURE() << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for
  } // for

  CINCH_CAPTURE() << separator;
  for(auto v: mesh_.vertices()) {
    auto xv = v->coordinates();
    CINCH_CAPTURE() << "^^^^ vertex id: " << v.id()
      << " with coordinates " << xv << endl;
    for(auto w: mesh_.wedges(v)) {

      CINCH_CAPTURE() << "     ---- wedge id: " << w.id() << endl;

      // Why do we have to do .front()?
      auto c = mesh_.cells(w).front();
      CINCH_CAPTURE() << "          ++++ cell id: " << c.id()
        << " with centroid " << c->centroid() << endl;

      auto e = mesh_.edges(w).front();
      CINCH_CAPTURE() << "          ++++ edge id: " << e.id()
        << " with midpoint " << e->midpoint() << endl;

      CINCH_CAPTURE() << "          ++++ side_facet_normal: "
        << w->facet_normal_left() << endl;

      CINCH_CAPTURE() << "          ++++ cell_facet_normal: "
        << w->facet_normal_right() << endl;

    } // for

  } // for

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the accessors
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, accessors) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "accessors" << endl;
  CINCH_CAPTURE() << separator;

  register_state(mesh_, "pressure", cells, real_t, persistent);
  register_state(mesh_, "density", cells, real_t);
  register_state(mesh_, "total energy", cells, real_t, persistent);
  register_state(mesh_, "velocity", vertices, vector_t, persistent);
  register_state(mesh_, "H", edges, vector_t);
  register_state(mesh_, "point_not_persistent", vertices, point_t);
  register_state(mesh_, "point_is_persistent", vertices, point_t, persistent);

  struct data_t {
    double x, y;
  };  
  register_global_state(mesh_, "const", data_t);

  CINCH_CAPTURE() << "Accessing state with type real_t:" << endl;

  auto vr = access_type(mesh_, real_t);
  for(auto v: vr) {
    CINCH_CAPTURE() << "\t" << v.label() << " has type real_t" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing state with type data_t:" << endl;

  auto vd = access_type(mesh_, data_t);
  for(auto v: vd) {
    CINCH_CAPTURE() << "\t" << v.label() << " has type data_t" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing state with type real_t at cells:" << endl;

  auto va = access_type_if(mesh_, real_t, is_at(cells));
  for(auto v: va) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type real_t and is at cells" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing persistent state with type real_t at cells:"
    << endl;

  auto vp = access_type_if(mesh_, real_t, is_persistent_at(cells));

  for(auto v: vp) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type real_t and is persistent at cells" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing state with type vector_t at vertices:" << endl;

  auto vv = access_type_if(mesh_, vector_t, is_at(vertices));
  for(auto v: vv) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type vector_t and is at vertices" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE()
    << "Accessing persistent state with type vector_t at vertices:" << endl;

  auto vpv = access_type_if(mesh_, vector_t, is_persistent_at(vertices));
  for(auto v: vpv) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type vector_t and is persistent at vertices" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing state with type vector_t at edges:" << endl;

  auto ve = access_type_if(mesh_, vector_t, is_at(edges));
  for(auto v: ve) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type vector_t and is at edges" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing state with type point_t at vertices:" << endl;

  auto pv = access_type_if(mesh_, point_t, is_at(vertices));
  for(auto v: pv) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type point_t and is at vertices" << endl;
  } // for

  CINCH_CAPTURE() << endl;

  CINCH_CAPTURE() << "Accessing persistent state with type point_t at vertices:"
            << endl;

  auto ppv = access_type_if(mesh_, point_t, is_persistent_at(vertices));
  for(auto v: ppv) {
    CINCH_CAPTURE() << "\t" << v.label() <<
      " has type point_t and is persistent at vertices" << endl;
  } // for

  CINCH_CAPTURE() << endl;

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the state
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, state) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "state" << endl;
  CINCH_CAPTURE() << separator;

  register_state(mesh_, "pressure", cells, real_t, persistent);
  register_state(mesh_, "velocity", vertices, vector_t, persistent);
  register_state(mesh_, "H", edges, vector_t);
  register_state(mesh_, "cornerdata", corners, int32_t);
  register_state(mesh_, "wedgedata", wedges, bool);

  struct data_t {
    int x, y;
  };  
  register_global_state(mesh_, "const", data_t);


  auto p = access_state(mesh_, "pressure", real_t);
  auto velocity = access_state(mesh_, "velocity", vector_t);
  auto H = access_state(mesh_, "H", vector_t);
  auto cd = access_state(mesh_, "cornerdata", int32_t);
  auto wd = access_state(mesh_, "wedgedata", bool);
  auto c = access_global_state(mesh_, "const", data_t);

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

  // corners
  ASSERT_EQ(16, mesh_.num_corners());
  for (auto c: mesh_.corners()) {
    cd[c] = c.id();
  } // for

  for (auto c: mesh_.corners()) {
    ASSERT_EQ(c.id(), cd[c]);
  } // for

  ASSERT_EQ(32, mesh_.num_wedges());
  // wedges
  for (auto w: wd) {
    wd[w] = (w%2) ? true : false;
  } // for

  for (auto w: wd) {
    if (w%2) {
      ASSERT_TRUE(wd[w]);
    }
    else {
      ASSERT_FALSE(wd[w]);
    }
  } // for

  // test global data
  c = { 1, 2 };
  ASSERT_EQ(1, c->x);  ASSERT_EQ(1, (*c).x);
  ASSERT_EQ(2, c->y);  ASSERT_EQ(2, (*c).y);

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief A final test to compare the blessed file and do CINCH_DUMP().
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, cinch_dump) {
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED("burton.blessed"));
  cout << CINCH_DUMP() << endl;
}

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
