/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// 
/// \brief Tests general features of the burton mesh.
////////////////////////////////////////////////////////////////////////////////


// user includes
#include "burton_3d_test.h"

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

  mesh_.dump( CINCH_CAPTURE() );

  cout << CINCH_DUMP() << endl;
  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));
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

  cout << CINCH_DUMP() << endl;
  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

} // TEST_F


////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh geometry functions
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_3d, geometry) {
  std::string separator;
  separator.insert(0,80,'=');
  separator.append("\n");

  CINCH_CAPTURE() << separator;
  CINCH_CAPTURE() << "geometry" << endl;

  CINCH_CAPTURE() << separator;
  for(auto c: mesh_.cells()) {
    auto xc = c->centroid();
    auto vol = c->volume();

    CINCH_CAPTURE() << "---- cell id: " << c.id()
      << " with centroid " << xc << " and volume " << vol << endl;

    for(auto v : mesh_.vertices(c)){
      auto xv = v->coordinates();
      CINCH_CAPTURE() << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for

  CINCH_CAPTURE() << separator;
  for(auto f: mesh_.faces()) {
    auto xc = f->centroid();
    auto a = f->area();
    auto n = f->normal();

    CINCH_CAPTURE() << "---- face id: " << f.id()
      << " with midpoint " << xc << ", area " << a
      << " and normal " << n << endl;

    for(auto v : mesh_.vertices(f)){
      auto xv = v->coordinates();
      CINCH_CAPTURE() << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for

  CINCH_CAPTURE() << separator;
  for(auto e: mesh_.edges()) {
    auto xc = e->midpoint();
    auto l = e->length();

    CINCH_CAPTURE() << "---- edge id: " << e.id()
      << " with midpoint " << xc << ", length " << l << endl;

    for(auto v : mesh_.vertices(e)){
      auto xv = v->coordinates();
      CINCH_CAPTURE() << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for

  cout << CINCH_DUMP() << endl;
  std::string outfile = output_prefix()+".blessed";
  CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));

} // TEST_F

////////////////////////////////////////////////////////////////////////////////
//! \brief test the face normals
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_3d, normals) {

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
    ASSERT_EQ( 3, mesh_.faces( cn ).size() );
    ASSERT_EQ( 3, mesh_.edges( cn ).size() );
    ASSERT_EQ( 1, mesh_.vertices( cn ).size() );
    ASSERT_EQ( 6, ws.size() );

    for ( auto wg = ws.begin(); wg != ws.end(); ++wg ) {
      // first edge
      {
        ASSERT_EQ( 1, mesh_.edges( *wg ).size() );
        ASSERT_EQ( 1, mesh_.faces( *wg ).size() );
        ASSERT_EQ( 1, mesh_.vertices( *wg ).size() );
        ASSERT_EQ( 1, mesh_.corners( *wg ).size() );
        auto fc = mesh_.faces(*wg).front();
        auto n = wg->facet_normal_right();
        auto fx = fc->centroid();
        auto cx = cl->centroid();
        auto delta = fx - cx;
        auto dot = dot_product( n, delta );
        ASSERT_GT( dot, 0 );
      }
      // second edge
      ++wg;
      {
        ASSERT_EQ( 1, mesh_.edges( *wg ).size() );
        ASSERT_EQ( 1, mesh_.faces( *wg ).size() );
        ASSERT_EQ( 1, mesh_.vertices( *wg ).size() );
        ASSERT_EQ( 1, mesh_.corners( *wg ).size() );
        auto fc = mesh_.faces(*wg).front();
        auto n = wg->facet_normal_left();
        auto fx = fc->centroid();
        auto cx = cl->centroid();
        auto delta = fx - cx;
        auto dot = dot_product( n, delta );
        ASSERT_GT( dot, 0 );
      }
    } // wedges

  } // for

  auto bnd_vertices = filter_boundary( mesh_.vertices() );

  for(auto vt : bnd_vertices) {
    
    auto ws = filter_boundary( mesh_.wedges(vt) );

    ASSERT_EQ( 0, ws.size()%2 );

    for ( auto wg = ws.begin(); wg != ws.end(); ++wg ) {
      // first edge
      {
        auto cl = mesh_.cells(*wg).front();
        auto fc = mesh_.faces(*wg).front();
        auto n = wg->facet_normal_right();
        auto fx = fc->centroid();
        auto cx = cl->centroid();
        auto delta = fx - cx;
        auto dot = dot_product( n, delta );
        ASSERT_GT( dot, 0 );
      }
      // second edge
      ++wg;
      {
        auto cl = mesh_.cells(*wg).front();
        auto fc = mesh_.faces(*wg).front();
        auto n = wg->facet_normal_left();
        auto fx = fc->centroid();
        auto cx = cl->centroid();
        auto delta = fx - cx;
        auto dot = dot_product( n, delta );
        ASSERT_GT( dot, 0 );
      }
    } // wedges

  }

} // TEST_F
