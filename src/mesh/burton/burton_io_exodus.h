/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  // 
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
/*!
 * \file io_exodus.h
 * \authors wohlbier
 * \date Initial file creation: Oct 07, 2015
 */

#pragma once

//! system includes
#include <cstring>

#ifdef HAVE_EXODUS
#  include <exodusII.h>
#endif

//! user includes
#include "flecsi/io/io_exodus.h"
#include "ale/mesh/burton/burton_mesh.h"


namespace flecsi {

//! bring burton mesh to front of namespace
using ale::mesh::burton_mesh_t;

//! Register file extension "g" with factory.
static bool burton_exodus_g_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType("g",
    flecsi::create_io_exodus<burton_mesh_t>);

//! Register file extension "exo" with factory.
static bool burton_exodus_exo_registered =
  flecsi::io_factory_t<burton_mesh_t>::instance().registerType("exo",
    flecsi::create_io_exodus<burton_mesh_t>);

/*!
  Implementation of exodus mesh read for burton specialization.

  \param[in] name Read burton mesh \e m from \e name.
  \param[out] m Populate burton mesh \e m with contents of \e name.

  \return Exodus error code. 0 on success.
 */
template<>
inline 
int32_t 
flecsi::io_exodus_t<burton_mesh_t>::read( const std::string &name,
                                          burton_mesh_t &m) 
{

#ifdef HAVE_EXODUS

  std::cout << "Reading mesh from: " << name << std::endl;

  // alias some types
  using std::vector;
  using   mesh_t = burton_mesh_t;
  using   real_t = typename mesh_t::real_t;
  using  point_t = typename mesh_t::point_t;
  using vector_t = typename mesh_t::vector_t;
  using vertex_t = typename mesh_t::vertex_t;
  using ex_real_t = real_t;
  using ex_index_t = int;

  // size of floating point variables used in app.
  int CPU_word_size = sizeof(real_t);
  // size of floating point stored in name.
  int IO_word_size = 0;
  float version;

  auto exoid = ex_open(
    name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version);
  assert(exoid >= 0);

  // get the initialization parameters
  ex_init_params exopar;
  auto status = ex_get_init_ext(exoid, &exopar);
  assert(status == 0);

  // verify mesh dimension
  assert(m.dimension() == exopar.num_dim);
  auto num_nodes = exopar.num_nodes;
  auto num_elem = exopar.num_elem;
  auto num_elem_blk = exopar.num_elem_blk;
  auto num_node_sets = exopar.num_node_sets;
  auto num_side_sets = exopar.num_side_sets;

  m.init_parameters(num_nodes);

  // read nodes
  vector<ex_real_t> xcoord(num_nodes);
  vector<ex_real_t> ycoord(num_nodes);
  status = ex_get_coord(exoid, xcoord.data(), ycoord.data(), nullptr);
  assert(status == 0);

  // put nodes into mesh
  vector<vertex_t *> vs;
  for (size_t i = 0; i < num_nodes; ++i) {
    point_t p = { static_cast<real_t>(xcoord[i]), 
                  static_cast<real_t>(ycoord[i]) };
    auto v = m.create_vertex(p);
    vs.push_back(v);
  } // for

  // 1 block for now

  // read blocks
  vector<ex_index_t> blockids( num_elem_blk );
  status = ex_get_elem_blk_ids(exoid, blockids.data() );
  assert(status == 0);

  char block_name[MAX_LINE_LENGTH];
  status = ex_get_name(exoid, EX_ELEM_BLOCK, blockids[0], block_name);
  assert(status == 0);

  // get the info about this block
  ex_index_t num_elem_this_blk = 0;
  ex_index_t num_attr = 0;
  ex_index_t num_nodes_per_elem = 0;
  char elem_type[MAX_STR_LENGTH];
  status = ex_get_elem_block(
    exoid, blockids[0], elem_type, &num_elem_this_blk, &num_nodes_per_elem, &num_attr);
  assert(status == 0);

  // verify that there is only one blk
  assert( num_elem_this_blk == num_elem );
  assert( num_elem_blk == 1 );

  // verify mesh has quads
  assert(num_nodes_per_elem == 4);
  assert(strcmp(elem_type, "quad") == 0);

  // read element definitions
  vector<ex_index_t> elt_conn(num_elem * num_nodes_per_elem);
  status = ex_get_elem_conn(exoid, blockids[0], elt_conn.data());
  assert(status == 0);

  // create cells in mesh
  for (size_t e = 0; e < num_elem; ++e) {
    auto b = e*num_nodes_per_elem; // base offset into elt_conn
    auto c = m.create_cell({vs[elt_conn[b+0]-1],
      vs[elt_conn[b+1]-1],
      vs[elt_conn[b+2]-1],
      vs[elt_conn[b+3]-1]});
  }
  m.init();

  // creating fields from an exodus file will be problematic since the type
  // information on the file has been thrown away. For example, an int field
  // is written to the file as a floating point field, and there is no way to
  // recover the fact that it is an int field. Furthermore, vector fields
  // are stored as individual scalar fields and recovering the fact that a
  // a collection of scalar fields makes a vector field would require clunky
  // processing of the variable names.


  return status;

#else

  std::cerr << "FLECSI not build with exodus support." << std::endl;
  std::exit(1);

  return -1;

#endif
}

/*!
  Implementation of exodus mesh write for burton specialization.

  \param[in] name Write burton mesh \e m to \e name.
  \param[in] m Burton mesh to write to \e name.

  \return Exodus error code. 0 on success.
 */
//FIXME: should allow for const mesh_t &
//int32_t io_exodus_t::write(
//    const std::string &name, const mesh_t &m) {
template<>
inline
int32_t 
flecsi::io_exodus_t<burton_mesh_t>::write( const std::string &name,
                                           burton_mesh_t &m) 
{

#ifdef HAVE_EXODUS

  std::cout << "Writing mesh to: " << name << std::endl;

  // alias some types
  using std::vector;
  using std::array;
  using   mesh_t = burton_mesh_t;
  using   real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t;
  using ex_real_t = real_t;
  using ex_index_t = int;


  //============================================================================
  // initial setup
  //============================================================================

  // size of floating point variables used in app.
  int CPU_word_size = sizeof(real_t);
  // size of floating point to be stored in file.
  int IO_word_size = sizeof(ex_real_t);
  // determine the file creation mode
  int cmode = EX_CLOBBER;

  auto exoid =
    ex_create(name.c_str(), cmode, &CPU_word_size, &IO_word_size);
  assert(exoid >= 0);

  // get the general statistics
  auto d = m.dimension();
  auto num_nodes = m.num_vertices();
  auto num_elem = m.num_cells();
  auto num_elem_blk = 1;
  auto num_node_sets = 0;
  auto num_side_sets = 0;

  // initialize the file.
  auto status = ex_put_init(exoid, "Exodus II output from flecsi.", d,
    num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets);
  assert(status == 0);

  // get the coordinates from the mesh.
  vector<ex_real_t> xcoord(num_nodes);
  vector<ex_real_t> ycoord(num_nodes);
  for (auto v : m.vertices()) {
    xcoord[v.id()] = v->coordinates()[0];
    ycoord[v.id()] = v->coordinates()[1];
  } // for
  // write the coordinates to the file
  status = ex_put_coord(exoid, xcoord.data(), ycoord.data(), nullptr);
  assert(status == 0);

  // write the coordinate names
  const char *coord_names[3];
  coord_names[0] = "x";
  coord_names[1] = "y";
  coord_names[2] = "z";
  status = ex_put_coord_names(exoid, (char **)coord_names);
  assert(status == 0);

  // loop over element blocks
  auto blockid = 0;
  auto num_attr = 0;
  auto num_nodes_per_elem = 4;
  status = ex_put_elem_block(
    exoid, blockid, "quad", num_elem, num_nodes_per_elem, num_attr);
  assert(status == 0);

  // element definitions
  vector<ex_index_t> elt_conn(num_elem * num_nodes_per_elem);
  auto i = 0;
  for (auto c : m.cells()) {
    for (auto v : m.vertices(c)) {
      elt_conn[i] = v.id() + 1; // 1 based index in exodus
      i++;
    } // for
  } // for

  // write connectivity
  status = ex_put_elem_conn(exoid, blockid, elt_conn.data());
  assert(status == 0);

  //============================================================================
  // write field data
  //============================================================================

  //----------------------------------------------------------------------------
  // nodal field data
  //----------------------------------------------------------------------------

  int num_nf = 0; // number of nodal fields
  // real scalars persistent at vertices
  auto rspav = access_type_if(m, real_t, is_persistent_at(vertices));
  num_nf += rspav.size();
  // int scalars persistent at vertices
  auto ispav = access_type_if(m, int, is_persistent_at(vertices));
  num_nf += ispav.size();
  // real vectors persistent at vertices
  auto rvpav = access_type_if(m, vector_t, is_persistent_at(vertices));
  num_nf += m.dimension()*rvpav.size();

  // variable extension for vectors
  std::string var_ext[3];
  var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

  // put the number of nodal fields
  if (num_nf > 0) {
    status = ex_put_var_param(exoid, "n", num_nf);
    assert(status == 0);
  }

  // fill node variable names array
  int inum = 1;
  for(auto sf: rspav) {
    status = ex_put_var_name(exoid, "n", inum++, sf.label().c_str());
    assert(status == 0);
  } // for
  for(auto sf: ispav) {
    status = ex_put_var_name(exoid, "n", inum++, sf.label().c_str());
  } // for
  for(auto vf: rvpav) {
    auto label = vf.label();
    for(int d=0; d < m.dimension(); ++d) {
      auto dim_label = label + var_ext[d];
      status = ex_put_var_name(exoid, "n", inum++, dim_label.c_str());
    } // for
  } // for


  // write nodal fields
  inum = 1;
  // node field buffer
  vector<ex_real_t> nf(num_nodes);
  for(auto sf: rspav) {
    for(auto v: m.vertices()) nf[v.id()] = sf[v];
    status = ex_put_nodal_var(exoid, /* timestep */ 1, inum++, num_nodes, nf.data());
    assert(status == 0);
  } // for
  for(auto sf: ispav) {
    // cast int fields to real_t
    for(auto v: m.vertices()) nf[v.id()] = (real_t)sf[v];
    status = ex_put_nodal_var(exoid, /* timestep */ 1, inum++, num_nodes, nf.data());
    assert(status == 0);
  } // for
  for(auto vf: rvpav) {
    for(int d=0; d < m.dimension(); ++d) {
      for(auto v: m.vertices()) nf[v.id()] = vf[v][d];
      status = ex_put_nodal_var(exoid, /* timestep */ 1, inum++, num_nodes, nf.data());
      assert(status == 0);
    } // for
  } // for



  //----------------------------------------------------------------------------
  // element field data
  //----------------------------------------------------------------------------

  int num_ef = 0; // number of element fields
  // real scalars persistent at cells
  auto rspac = access_type_if(m, real_t, is_persistent_at(cells));
  num_ef += rspac.size();
  // int scalars persistent at cells
  auto ispac = access_type_if(m, int, is_persistent_at(cells));
  num_ef += ispac.size();
  // real vectors persistent at cells
  auto rvpac = access_type_if(m, vector_t, is_persistent_at(cells));
  num_ef += m.dimension()*rvpac.size();


  // put the number of element fields
  if(num_ef > 0) {
    status = ex_put_var_param(exoid, "e", num_ef);
    assert(status == 0);
  } // if

  // fill element variable names array
  inum = 1;
  for(auto sf: rspac) {
    status = ex_put_var_name(exoid, "e", inum++, sf.label().c_str());
    assert(status == 0);
  } // for
  for(auto sf: ispac) {
    status = ex_put_var_name(exoid, "e", inum++, sf.label().c_str());
    assert(status == 0);
  } // for
  for(auto vf: rvpac) {
    auto label = vf.label();
    for(int d=0; d < m.dimension(); ++d) {
      auto dim_label = label + var_ext[d];
      status = ex_put_var_name(exoid, "e", inum++, dim_label.c_str());
    } // for
  } // for

  // write element fields
  inum = 1;
  // element field buffer
  vector<ex_real_t> ef(num_elem);
  for(auto sf: rspac) {
    for(auto c: m.cells()) ef[c.id()] = sf[c];
    status = ex_put_elem_var(exoid, /* timestep */ 1, inum++, /* blkid */ 0, num_elem, ef.data());
    assert(status == 0);
  } // for
  for(auto sf: ispac) {
    // cast int fields to real_t
    for(auto c: m.cells()) ef[c.id()] = (real_t)sf[c];
    status = ex_put_elem_var(exoid, /* timestep */ 1, inum++, /* blkid */ 0, num_elem, ef.data());
    assert(status == 0);
  } // for
  for(auto vf: rvpac) {
    for(int d=0; d < m.dimension(); ++d) {
      for(auto c: m.cells()) ef[c.id()] = vf[c][d];
      status = ex_put_elem_var(exoid, /* timestep */ 1, inum++, /* blkid */ 0, num_elem, ef.data());
      assert(status == 0);
    } // for
  } // for


  //============================================================================
  // final setup
  //============================================================================
  // set the time
  ex_real_t soln_time = m.get_time();
  status = ex_put_time(exoid, /* timestep */ 1, &soln_time );
  assert(status == 0);

  // close
  status = ex_close(exoid);

  return status;


#else

  std::cerr << "FLECSI not build with exodus support." << std::endl;
  std::exit(1);

  return -1;

#endif

} // io_exodus_t::write

} // namespace flecsi

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
