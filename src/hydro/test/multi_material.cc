

// create some field data
register_state( mesh, "density",   cells,   real_t, persistent );
register_state( mesh, "pressure",  cells,   real_t, persistent );
register_state( mesh, "velocity", vertex, vector_t, persistent );

// register materials
register_materials( mesh, { "iron", "copper" } );

// register material properties.  these are sparse!
register_material_state( mesh, "mass",      cells,   real_t, persistent );
register_material_state( mesh, "volume",    cells,   real_t, persistent );
register_material_state( mesh, "density",   cells,   real_t, persistent );
register_material_state( mesh, "temperature", cells,   real_t, persistent );
register_material_state( mesh, "velocity", vertex, vector_t, persistent );

// get the list of materials
auto mats = registered_materials( mesh );
auto num_mats = mats.size();

// get global field accessors
auto d = access_state("density", real_t);
auto u = access_state("velocity", vector_t);

// get material field accessors
auto mi = access_material_state("mass",       real_t);
auto vi = access_material_state("volume",     real_t);
auto di = access_material_state("density",    real_t);
auto ti = access_material_state("temperature", real_t);
auto ui = access_material_state("velocity", vector_t);

//-----------------------------------------------------------------------------
// something cell-based
for ( auto cell : mesh.cells() ) {
  d[cell] = 0; // access a global cell field

  // the materials in this cell
  auto cell_mats = mesh.materials(cell);

  // gather some material info to a cell quantity
  for ( auto mat : cell_mats ) 
    d[cell] += mi[mat] / vi[mat]; // mat is really a 2d index space ->
                                  // cell id + material id

  d[cell] /= cell_mats.size();  
 }


//-----------------------------------------------------------------------------
// something vertex based
for ( auto vert : mesh.vertices() ) {
  u[vert] = 0; // access a vertex based field
  
  // the materials assigned to this vertex.  this is different
  // than the assigned cell materials!!!
  auto vert_mats = mesh.materials(vert);

  // gather some material info to a vertex quantity
  for ( auto mat : vert_mats ) 
    u[vert] += ui[mat]; // remember, mat is a 2d index space
  
  d[c] /= mats.size();  
 }


//-----------------------------------------------------------------------------
// something over the whole sparse storage

for ( auto i : di.index_space() )
  di[i] = mi[i] / vi[i];

//-----------------------------------------------------------------------------
// something material based
for ( auto mat : mats ) {

  // get this materials eos
  auto eos = eos_list[ mat->id() ];

  // now apply it to something
  for ( auto cell : mesh.cells(mat) )
    ti[cell] = eos.get_temperature( di[cell] );  // cell is now really
                                                 // a 2d index space
    
 }


