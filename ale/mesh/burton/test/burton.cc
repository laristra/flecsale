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

  for(auto v : mesh_.vertices()) {
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Edges in mesh:" << endl;

  for(auto e : mesh_.edges()) {
    CINCH_CAPTURE() << "----------- edge id: " << e.id()
      << " with midpoint " << e->midpoint() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Cells in mesh:" << endl;

  for(auto c : mesh_.cells()) {
    CINCH_CAPTURE() << "----------- cell id: " << c.id()
      << " with centroid " << c->centroid() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Corners in mesh:" << endl;

  for(auto c : mesh_.corners()) {
    CINCH_CAPTURE() << "----------- corner id: " << c.id() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "Wedges in mesh:" << endl;

  for(auto w : mesh_.wedges()) {
    CINCH_CAPTURE() << "----------- wedge id: " << w.id() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each vertex:" << endl;

  for(auto v: mesh_.vertices()) {
    CINCH_CAPTURE() << "^^^^^^^^Vertex id: " << v.id() << endl;

    CINCH_CAPTURE() << "    ----Corners:" << endl;
    for(auto c: mesh_.corners(v))
      CINCH_CAPTURE() << "    ++++ corner id: " << c.id() << endl;

    CINCH_CAPTURE() << "    ----Cells:" << endl;
    for(auto c: mesh_.cells(v))
      CINCH_CAPTURE() << "    ++++ cell id: " << c.id() << endl;

    CINCH_CAPTURE() << "    ----Edges:" << endl;
    for(auto e: mesh_.edges(v))
      CINCH_CAPTURE() << "    ++++ edge id: " << e.id() << endl;
  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each edge:" << endl;

  for(auto e : mesh_.edges()) {
    CINCH_CAPTURE() << "^^^^^^^^Edge id: " << e.id() << endl;

    CINCH_CAPTURE() << "    ----Corners:" << endl;
    for(auto cnr : mesh_.corners(e))
      CINCH_CAPTURE() << "    ++++ corner id: " << cnr.id() << endl;

    CINCH_CAPTURE() << "    ----Cells:" << endl;
    for(auto c : mesh_.cells(e))
      CINCH_CAPTURE() << "    ++++ cell id: " << c.id() << endl;

    CINCH_CAPTURE() << "    ----Vertices:" << endl;
    for(auto v : mesh_.vertices(e))
      CINCH_CAPTURE() << "    ++++ vertex id: " << v.id() << endl;

  } // for

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each cell:" << endl;

  for(auto c : mesh_.cells()) {
    CINCH_CAPTURE() << "^^^^^^^^Cell id: " << c.id() << endl;

    CINCH_CAPTURE() << "    ----Corners:" << endl;
    for(auto cnr : mesh_.corners(c))
      CINCH_CAPTURE() << "    ++++ corner id: " << cnr.id() << endl;

    CINCH_CAPTURE() << "    ----Edges:" << endl;
    for(auto e : mesh_.edges(c))
      CINCH_CAPTURE() << "    ++++ edge id: " << e.id() << endl;

    CINCH_CAPTURE() << "    ----Vertices:" << endl;
    for(auto v : mesh_.vertices(c))
      CINCH_CAPTURE() << "    ++++ vertex id: " << v.id() << endl;

  } // for


  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each corner:" << endl;

  for(auto c : mesh_.corners()) {
    CINCH_CAPTURE() << "^^^^^^^^Corner id: " << c.id() << endl;

#if 0
    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w: mesh_.wedges(c)) 
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;
#endif

    CINCH_CAPTURE() << "    ----Cells:" << endl;
    for(auto cl: mesh_.cells(c)) 
      CINCH_CAPTURE() << "    ++++ cell id: " << cl.id() << endl;

    CINCH_CAPTURE() << "    ----Edges:" << endl;
    for(auto e: mesh_.edges(c)) 
      CINCH_CAPTURE() << "    ++++ edge id: " << e.id() << endl;

    CINCH_CAPTURE() << "    ----Vertices:" << endl;
    for(auto v: mesh_.vertices(c)) 
      CINCH_CAPTURE() << "    ++++ vertex id: " << v.id() << endl;

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
    auto area = c->area();

    CINCH_CAPTURE() << "---- cell id: " << c.id()
      << " with centroid " << xc << " and area " << area << endl;

    for(auto v : mesh_.vertices(c)){
      auto xv = v->coordinates();
      CINCH_CAPTURE() << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for


  CINCH_CAPTURE() << separator;
  for(auto e: mesh_.edges()) {
    auto xc = e->midpoint();
    auto l = e->length();
    auto n = e->normal();

    CINCH_CAPTURE() << "---- edge id: " << e.id()
      << " with midpoint " << xc << ", length " << l 
      << " and normal " << n << endl;

    for(auto v : mesh_.vertices(e)){
      auto xv = v->coordinates();
      CINCH_CAPTURE() << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for


  CINCH_CAPTURE() << separator;
  for(auto c: mesh_.cells()) {
    CINCH_CAPTURE() << "^^^^^^^^Cell id: " << c.id() << endl;

    vector_t sum(0);

    CINCH_CAPTURE() << "    ----Corners:" << endl;
    for ( auto cn: mesh_.corners(c) ) {
      
     CINCH_CAPTURE() << "    ++++ corner id: " << cn.id() << endl;

      for ( auto w: cn->wedges() ) {
        auto left  = w->facet_normal_left();
        auto right = w->facet_normal_right();
        CINCH_CAPTURE() << "    .... wedge" << endl;
        CINCH_CAPTURE() << "         facet_normal_left: "  << left << endl;
        CINCH_CAPTURE() << "         facet_normal_right: " << right << endl;
        sum += left;
        sum += right;
      } // for

    } // for

    ASSERT_NEAR( 0, sum[0], test_tolerance );
    ASSERT_NEAR( 0, sum[1], test_tolerance );
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

  CINCH_CAPTURE() << "Accessing state with type real_t" << endl;

  std::vector<std::string> labels;
  auto vr = access_type(mesh_, real_t);
  for(auto v: vr) {
    std::cout << v.label() << std::endl;
    labels.push_back(v.label());
  } // for

  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "pressure")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "density")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "total energy")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing state with type data_t" << endl;

  auto vd = access_type(mesh_, data_t);
  for(auto v: vd) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "const")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing state with type real_t at cells" << endl;

  auto va = access_type_if(mesh_, real_t, is_at(cells));
  for(auto v: va) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "pressure")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "density")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "total energy")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing persistent state with type real_t at cells"
    << endl;

  auto vp = access_type_if(mesh_, real_t, is_persistent_at(cells));

  for(auto v: vp) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "pressure")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "total energy")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing state with type vector_t at vertices" << endl;

  auto vv = access_type_if(mesh_, vector_t, is_at(vertices));
  for(auto v: vv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "velocity")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE()
    << "Accessing persistent state with type vector_t at vertices" << endl;

  auto vpv = access_type_if(mesh_, vector_t, is_persistent_at(vertices));
  for(auto v: vpv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "velocity")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing state with type vector_t at edges" << endl;

  auto ve = access_type_if(mesh_, vector_t, is_at(edges));
  for(auto v: ve) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "H")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing state with type point_t at vertices" << endl;

  auto pv = access_type_if(mesh_, point_t, is_at(vertices));
  for(auto v: pv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "point_is_persistent")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "point_not_persistent")
    != labels.end());
  labels.clear();

  CINCH_CAPTURE() << "Accessing persistent state with type point_t at vertices"
            << endl;

  auto ppv = access_type_if(mesh_, point_t, is_persistent_at(vertices));
  for(auto v: ppv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "point_is_persistent")
    != labels.end());
  labels.clear();

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
  register_state(mesh_, "cornerdata", corners, integer_t);
  register_state(mesh_, "wedgedata", wedges, bool);

  struct data_t {
    int x, y;
  };  
  register_global_state(mesh_, "const", data_t);


  auto p = access_state(mesh_, "pressure", real_t);
  auto velocity = access_state(mesh_, "velocity", vector_t);
  auto H = access_state(mesh_, "H", vector_t);
  auto cd = access_state(mesh_, "cornerdata", integer_t);
  auto wd = access_state(mesh_, "wedgedata", bool);

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
  auto cnst = access_global_state(mesh_, "const", data_t);
  cnst = { 1, 2 };
  ASSERT_EQ(1, cnst->x);  ASSERT_EQ(1, (*cnst).x);
  ASSERT_EQ(2, cnst->y);  ASSERT_EQ(2, (*cnst).y);

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief A final test to compare the blessed file and do CINCH_DUMP().
////////////////////////////////////////////////////////////////////////////////
TEST_F(Burton, cinch_dump) {
  cout << CINCH_DUMP() << endl;
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED("burton.blessed"));
}

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
