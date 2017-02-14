/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some functionality for creating meshes.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#ifdef HAVE_OPENSSL
#  include <flecsi/utils/checksum.h>
#endif

// system includes
#include <iomanip>

namespace ale {
namespace mesh {

////////////////////////////////////////////////////////////////////////////////
//! \brief Output the checksums of solution quantities.
//!
//! \param [in] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t checksum( T & mesh ) 
{

#ifdef HAVE_OPENSSL

  using mesh_t = T;
  using counter_t = typename mesh_t::counter_t;
  using integer_t = typename mesh_t::integer_t;
  using real_t = typename mesh_t::real_t;
  using vector_t = typename mesh_t::vector_t; 

  constexpr auto num_dims = mesh_t::num_dimensions;

  //----------------------------------------------------------------------------
  // create some lambda functions
  auto checksum = []( auto && vs, auto && f )
  {
    using field_type = std::decay_t< decltype( f[0] ) >;
    auto num_verts = vs.size();
    std::vector< field_type > vals(num_verts);
    #pragma omp parallel for
    for ( counter_t i=0; i<num_verts; ++i ) vals[i] = f[vs[i]];
    flecsi::utils::checksum_t cs;
    flecsi::utils::checksum(vals.data(), num_verts, cs);
    return cs;
  };

  // variable extension for vectors
  std::string var_ext[3];
  var_ext[0] = "_x"; var_ext[1] = "_y";  var_ext[2] = "_z";

  std::cout << std::string(65, '=') << std::endl;
  std::cout << std::left << std::setw(32) << "Field Name"  << " " 
            << std::setw(32) << std::left << "Checksum" << std::endl;
  std::cout << std::string(65, '-') << std::endl;

  //----------------------------------------------------------------------------
  // Checksum for Coordinates
  auto verts = mesh.vertices();
  auto num_verts = verts.size();

  // get the coordinates from the mesh.
  for(int d=0; d < num_dims; ++d) {
    std::vector< real_t > vals(num_verts);
    #pragma omp parallel for
    for ( counter_t i=0; i<num_verts; ++i ) vals[i] = verts[i]->coordinates()[d];
    flecsi::utils::checksum_t cs;
    flecsi::utils::checksum(vals.data(), num_verts, cs);
    std::cout << std::left << std::setw(32) << "node_coordinates"+var_ext[d] 
              << " " << std::setw(32) << cs.strvalue << std::endl;

  } // for

  //----------------------------------------------------------------------------
  // Checksum Nodal Solution Quantities

  // real scalars persistent at vertices
  auto rspav = get_accessors_all(
    mesh, real_t, dense, 0, has_attribute_at(persistent,vertices)
  );
  for(auto sf: rspav) {
    auto cs = checksum( verts, sf );
    std::cout << std::left << std::setw(32) << sf.label() << " " 
              << std::setw(32) << cs.strvalue << std::endl;
  }

  // int scalars persistent at vertices
  auto ispav = get_accessors_all(
    mesh, integer_t, dense, 0, has_attribute_at(persistent,vertices)
  );
  for(auto sf: ispav) {
    auto cs = checksum( verts, sf );
    std::cout << std::left << std::setw(32) << sf.label() << " " 
              << std::setw(32) << cs.strvalue << std::endl;
  }

  // real vectors persistent at vertices
  auto rvpav = get_accessors_all(
    mesh, vector_t, dense, 0, has_attribute_at(persistent,vertices)
  );
  for(auto vf: rvpav) {
    for(int d=0; d < num_dims; ++d) {
      auto cs = checksum( verts, vf );
      std::cout << std::left << std::setw(32) << vf.label()+var_ext[d] 
                << " " << std::setw(32) << cs.strvalue << std::endl;
    } // for
  } // for


  //----------------------------------------------------------------------------
  // Checksum Cell Solution Quantities
  auto cels = mesh.cells();

  // real scalars persistent at cells
  auto rspac = get_accessors_all(
    mesh, real_t, dense, 0, has_attribute_at(persistent,cells)
  );
  for(auto sf: rspac) {
    auto cs = checksum( cels, sf );
    std::cout << std::left << std::setw(32) << sf.label() << " " 
              << std::setw(32) << cs.strvalue << std::endl;
  }
  // int scalars persistent at cells
  auto ispac = get_accessors_all(
    mesh, integer_t, dense, 0, has_attribute_at(persistent,cells)
  );
  for(auto sf: ispac) {
    auto cs = checksum( cels, sf );
    std::cout << std::left << std::setw(32) << sf.label() << " " 
              << std::setw(32) << cs.strvalue << std::endl;
  }
  // real vectors persistent at cells
  auto rvpac = get_accessors_all(
    mesh, vector_t, dense, 0, has_attribute_at(persistent,cells)
  );
  for(auto vf: rvpac) {
    for(int d=0; d < num_dims; ++d) {
      auto cs = checksum( cels, vf );
      std::cout << std::left << std::setw(32) << vf.label()+var_ext[d] << " " 
                << std::setw(32) << cs.strvalue << std::endl;
    } // for
  } // for

  std::cout << std::string(65, '=') << std::endl;


#endif // ENABLE_OPENSSL

  return 0;
}

} // namespace
} // namespace
