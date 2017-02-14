/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
//
// \file
// 
// \brief Tests general features of the burton mesh.
//
////////////////////////////////////////////////////////////////////////////////

// user includes
#include "burton_2d_test.h"

// using statements
using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////////////
//! \brief dump the mesh to std out
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, dump) {

  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  mesh_.dump( CINCH_CAPTURE() );

  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));
}

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh connectivity
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, mesh) {

  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

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

    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w: mesh_.wedges(v))
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;

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

    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w: mesh_.wedges(e))
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;

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

    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w: mesh_.wedges(c))
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;

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

    CINCH_CAPTURE() << "    ----Wedges:" << endl;
    for(auto w: mesh_.wedges(c)) 
      CINCH_CAPTURE() << "    ++++ wedge id: " << w.id() << endl;

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

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "For each wedge:" << endl;

  for(auto w : mesh_.wedges()) {
    CINCH_CAPTURE() << "^^^^^^^^Wedge id: " << w.id() << endl;

    CINCH_CAPTURE() << "    ----Corners:" << endl;
    for(auto c: mesh_.corners(w)) 
      CINCH_CAPTURE() << "    ++++ corner id: " << c.id() << endl;

    CINCH_CAPTURE() << "    ----Cells:" << endl;
    for(auto cl: mesh_.cells(w)) 
      CINCH_CAPTURE() << "    ++++ cell id: " << cl.id() << endl;

    CINCH_CAPTURE() << "    ----Edges:" << endl;
    for(auto e: mesh_.edges(w)) 
      CINCH_CAPTURE() << "    ++++ edge id: " << e.id() << endl;

    CINCH_CAPTURE() << "    ----Vertices:" << endl;
    for(auto v: mesh_.vertices(w)) 
      CINCH_CAPTURE() << "    ++++ vertex id: " << v.id() << endl;

  } // for

  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh geometry functions
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, geometry) {

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

      for ( auto w: mesh_.wedges(cn) ) {
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

  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the face normals
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, normals) {

  for(auto f : mesh_.faces()) {
    auto n = f->normal();
    auto fx = f->centroid();
    auto c = mesh_.cells(f).front();
    auto cx = c->centroid();
    auto delta = fx - cx;
    auto dot = dot_product( n, delta );
    ASSERT_GT( dot, 0 );
  } // for

  for(auto cn : mesh_.corners()) {
    auto cl = mesh_.cells(cn).front();
    auto vt = mesh_.vertices(cn).front();
    auto ws = mesh_.wedges(cn);
    ASSERT_EQ( 1, mesh_.cells( cn ).size() );
    ASSERT_EQ( 2, mesh_.edges( cn ).size() );
    ASSERT_EQ( 1, mesh_.vertices( cn ).size() );
    ASSERT_EQ( 2, ws.size() );
    // first edge
    {
      ASSERT_EQ( 1, mesh_.edges( ws[0] ).size() );
      ASSERT_EQ( 1, mesh_.vertices( ws[0] ).size() );
      ASSERT_EQ( 1, mesh_.corners( ws[0] ).size() );
      auto e = mesh_.edges(ws[0]).front();
      auto n = ws[0]->facet_normal_right();
      auto ex = e->midpoint();
      auto cx = cl->centroid();
      auto delta = ex - cx;
      auto dot = dot_product( n, delta );
      ASSERT_GT( dot, 0 );
    }
    // second edge
    {
      ASSERT_EQ( 1, mesh_.edges( ws[1] ).size() );
      ASSERT_EQ( 1, mesh_.vertices( ws[1] ).size() );
      ASSERT_EQ( 1, mesh_.corners( ws[1] ).size() );
      auto e = mesh_.edges(ws[1]).front();
      auto n = ws[1]->facet_normal_left();
      auto ex = e->midpoint();
      auto cx = cl->centroid();
      auto delta = ex - cx;
      auto dot = dot_product( n, delta );
      ASSERT_GT( dot, 0 );
    }
  } // for

} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test the accessors
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, accessors) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  cout << separator;
  cout << "accessors" << endl;
  cout << separator;

  register_data(mesh_, hydro, pressure, real_t, dense, 1, cells);
  register_data(mesh_, hydro, density, real_t, dense, 1, cells);
  register_data(mesh_, hydro, total_energy, real_t, dense, 1, cells);
  register_data(mesh_, hydro, velocity, vector_t, dense, 1, vertices);
  register_data(mesh_, hydro, H, vector_t, dense, 1, edges);
  register_data(mesh_, hydro, point_not_persistent, point_t, dense, 1, vertices);
  register_data(mesh_, hydro, point_is_persistent, point_t, dense, 1, vertices);

  get_accessor(mesh_, hydro, pressure, real_t, dense, 0).attributes().set( persistent );
  get_accessor(mesh_, hydro, total_energy, real_t, dense, 0).attributes().set( persistent );
  get_accessor(mesh_, hydro, velocity, vector_t, dense, 0).attributes().set( persistent );
  get_accessor(mesh_, hydro, point_is_persistent, point_t, dense, 0).attributes().set( persistent );

  struct data_t {
    double x, y;
  };  
  register_data(mesh_, hydro, const, data_t, global, 1);

  cout << "Accessing state with type real_t" << endl;

  std::vector<std::string> labels;
  auto vr = get_accessors(mesh_, hydro, real_t, dense, 0);
  for(auto v: vr) {
    std::cout << v.label() << std::endl;
    labels.push_back(v.label());
  } // for

  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "pressure")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "density")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "total_energy")
    != labels.end());
  labels.clear();

  cout << "Accessing state with type data_t" << endl;

  auto vd = get_accessors(mesh_, hydro, data_t, global, 0);
  for(auto v: vd) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "const")
    != labels.end());
  labels.clear();

  cout << "Accessing state with type real_t at cells" << endl;

  auto va = get_accessors(mesh_, hydro, real_t, dense, 0, is_at(cells));
  for(auto v: va) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "pressure")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "density")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "total_energy")
    != labels.end());
  labels.clear();

  cout << "Accessing persistent state with type real_t at cells"
    << endl;

  auto vp = get_accessors(mesh_, hydro, real_t, dense, 0, has_attribute_at(persistent, cells));

  for(auto v: vp) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "pressure")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "total_energy")
    != labels.end());
  labels.clear();

  cout << "Accessing state with type vector_t at vertices" << endl;

  auto vv = get_accessors(mesh_, hydro, vector_t, dense, 0, is_at(vertices));
  for(auto v: vv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "velocity")
    != labels.end());
  labels.clear();

  cout
    << "Accessing persistent state with type vector_t at vertices" << endl;

  auto vpv = get_accessors(mesh_, hydro, vector_t, dense, 0, has_attribute_at(persistent, vertices));
  for(auto v: vpv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "velocity")
    != labels.end());
  labels.clear();

  cout << "Accessing state with type vector_t at edges" << endl;

  auto ve = get_accessors(mesh_, hydro, vector_t, dense, 0, is_at(edges));
  for(auto v: ve) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "H")
    != labels.end());
  labels.clear();

  cout << "Accessing state with type point_t at vertices" << endl;

  auto pv = get_accessors(mesh_, hydro, point_t, dense, 0, is_at(vertices));
  for(auto v: pv) {
    labels.push_back(v.label());
  } // for
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "point_is_persistent")
    != labels.end());
  ASSERT_TRUE(std::find(labels.begin(), labels.end(), "point_not_persistent")
    != labels.end());
  labels.clear();

  cout << "Accessing persistent state with type point_t at vertices"
            << endl;

  auto ppv = get_accessors(mesh_, hydro, point_t, dense, 0, has_attribute_at(persistent,vertices));
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
TEST_F(burton_2d, state) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  cout << separator;
  cout << "state" << endl;
  cout << separator;

  register_data(mesh_, hydro, pressure, real_t, dense, 1, cells);
  register_data(mesh_, hydro, velocity, vector_t, dense, 1, vertices);
  register_data(mesh_, hydro, H, vector_t, dense, 1, edges);
  register_data(mesh_, hydro, cornerdata, integer_t, dense, 1, corners);
  register_data(mesh_, hydro, wedgedata, bool, dense, 1, wedges);

  get_accessor(mesh_, hydro, pressure, real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh_, hydro, velocity, vector_t, dense, 0).attributes().set(persistent);

  struct data_t {
    int x, y;
  };  
  register_data(mesh_, hydro, const, data_t, global, 1);


  auto p = get_accessor(mesh_, hydro, pressure, real_t, dense, 0);
  auto velocity = get_accessor(mesh_, hydro, velocity, vector_t, dense, 0);
  auto H = get_accessor(mesh_, hydro, H, vector_t, dense, 0);
  auto cd = get_accessor(mesh_, hydro, cornerdata, integer_t, dense, 0);
  auto wd = get_accessor(mesh_, hydro, wedgedata, bool, dense, 0);

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
  auto cnst = get_accessor(mesh_, hydro, const, data_t, global, 0);
  *cnst = { 1, 2 };
  ASSERT_EQ(1, cnst->x);  ASSERT_EQ(1, (*cnst).x);
  ASSERT_EQ(2, cnst->y);  ASSERT_EQ(2, (*cnst).y);

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test copying the mesh
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, copy) {

  for(auto v : mesh_.vertices()){
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << std::endl;
  } // for

  auto mesh_copy = mesh_;

  for(auto v : mesh_copy.vertices()){
    CINCH_CAPTURE() << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << std::endl;
  } // for

  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

} // TEST
