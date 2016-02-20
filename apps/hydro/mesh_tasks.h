
#pragma once

namespace apps {
namespace hydro {

////////////////////////////////////////////////////////////////////////////////
//! \brief The main driver for initializing the mesh
//!
//! \param [in,out] mesh the mesh object
//! \return 0 for success
////////////////////////////////////////////////////////////////////////////////
template< typename T >
int32_t init_mesh(T & mesh) 
{


  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 2;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.2;

  constexpr real_t x0 = -length_x/2.0;
  constexpr real_t y0 = -length_y/2.0;


  // reserve storage for the mesh
  auto num_vertex = ( num_cells_x + 1 ) * ( num_cells_y + 1 );
  mesh.init_parameters( num_vertex );
  
  
  // create the individual vertices
  using vertex_t = mesh_t::vertex_t;
  vector<vertex_t*> vs;
  
  auto delta_x = length_x / num_cells_x;
  auto delta_y = length_y / num_cells_y;

  auto num_vert_x = num_cells_x + 1;
  auto num_vert_y = num_cells_y + 1;

  for(size_t j = 0; j < num_vert_y; ++j) {
    auto y = y0 + j*delta_y;
    for(size_t i = 0; i < num_vert_x; ++i) {
      auto x = x0 + i*delta_x;
      auto v = mesh.create_vertex( {x, y} );
      vs.push_back(v);
    }
    
  }
  
  // define each cell
  auto index = [=](auto i, auto j) { return i + num_vert_x*j; };
  
  for(size_t j = 0; j < num_cells_y; ++j)
    for(size_t i = 0; i < num_cells_x; ++i) {
      auto c = 
        mesh.create_cell({
              vs[ index(i  , j  ) ],
              vs[ index(i+1, j  ) ],
              vs[ index(i+1, j+1) ],
              vs[ index(i  , j+1) ]});
    }
  
  
  // now finalize the mesh setup
  mesh.init();

}

} // namespace hydro
} // namespace apps
