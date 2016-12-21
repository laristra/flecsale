/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro incdludes
#include "tasks.h"
#include "types.h"

// user includes
#include <ale/mesh/factory.h>
#include <ale/utils/time_utils.h>

// system includes
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

// everything is in the hydro namespace
using namespace apps::hydro;

// right now cases are hard coded
//#define SODX_2D
//#define SODX_3D
//#define SODY_2D
//#define SODXY_2D
//#define SEDOV_2D
#define SEDOV_3D

///////////////////////////////////////////////////////////////////////////////
//! \brief A sample test of the hydro solver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

  // enable exceptions
#if defined(_GNU_SOURCE) && !defined(NDEBUG)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif


  //===========================================================================
  // SOD-X Inputs
  //===========================================================================

#ifdef SODX_2D

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sodx_1d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[0] < 0.0 ) {
        d = 1.0;
        p = 1.0;
      }
      else {
        d = 0.125;
        p = 0.1;
      }    
      return std::make_tuple( d, v, p );
    };

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // create a single boundary right now
  auto bc_type = boundary_condition_t< mesh_t::num_dimensions >();
  
  // install each boundary
  // -  both +ve and -ve side boundaries can be installed at once since 
  //    they will never overlap

  // x-boundaries
  auto bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[0] == -0.5*length_x || fx[0] == 0.5*length_x );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // y-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[1] == -0.5*length_y || fx[1] == 0.5*length_y );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );

  
  //===========================================================================
  // SOD-Y Inputs
  //===========================================================================

#elif defined(SODY_2D)

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sody_2d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 10;
  constexpr size_t num_cells_y = 100;

  constexpr real_t length_x = 0.1;
  constexpr real_t length_y = 1.0;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[1] < 0.0 ) {
        d = 1.0;
        p = 1.0;
      }
      else {
        d = 0.125;
        p = 0.1;
      }    
      return std::make_tuple( d, v, p );
    };

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // create a single boundary right now
  auto bc_type = boundary_condition_t< mesh_t::num_dimensions >();
  
  // install each boundary
  // -  both +ve and -ve side boundaries can be installed at once since 
  //    they will never overlap

  // x-boundaries
  auto bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[0] == -0.5*length_x || fx[0] == 0.5*length_x );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // y-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[1] == -0.5*length_y || fx[1] == 0.5*length_y );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );

  //===========================================================================
  // SOD-XY Inputs
  //===========================================================================

#elif defined(SODXY_2D)

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sodxy_2d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;
  size_t max_steps = 1e6;

  // rotate by 45 degrees
  constexpr real_t degrees = 45;

  // this is a lambda function to set the initial conditions
  auto ics = [&degrees] ( const auto & x )
    {
      // compute some angles
      auto radians = degrees * math::pi / 180;
      auto cos = std::cos( radians );
      auto sin = std::sin( radians );

      // compute rotation matrix
      vector_t xrot;
      xrot[0] = x[0]*cos - x[1]*sin;
      xrot[1] = x[0]*sin + x[1]*cos;
      
      // compute solution
      real_t d, p;
      vector_t v(0);
      if ( xrot[1] < 0.0 ) {
        d = 1.0;
        p = 1.0;
      }
      else {
        d = 0.125;
        p = 0.1;
      }    
      return std::make_tuple( d, v, p );
    };

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( num_cells_x, num_cells_y, length_x, length_y );

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // create a single boundary right now
  auto bc_type = boundary_condition_t< mesh_t::num_dimensions >();
  
  // install each boundary
  // -  both +ve and -ve side boundaries can be installed at once since 
  //    they will never overlap

  // x-boundaries
  auto bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[0] == -0.5*length_x || fx[0] == 0.5*length_x );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // y-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[1] == -0.5*length_y || fx[1] == 0.5*length_y );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );

  // now rotate the mesh
  mesh::rotate( mesh, degrees );


  //===========================================================================
  // SOD-X Inputs
  //===========================================================================

#elif defined(SODX_3D)

  using mesh_t = mesh_3d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix = "sodx_3d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 1;

  // the grid dimensions
  constexpr size_t num_cells_x = 100;
  constexpr size_t num_cells_y = 10;
  constexpr size_t num_cells_z = 10;

  constexpr real_t length_x = 1.0;
  constexpr real_t length_y = 0.1;
  constexpr real_t length_z = 0.1;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 1.0, .volume = 1.0, .growth = 1.e3 };
  constexpr real_t final_time = 0.2;
  constexpr real_t initial_time_step = 1.0;
  size_t max_steps = 1e6;

  // this is a lambda function to set the initial conditions
  auto ics = [] ( const auto & x )
    {
      real_t d, p;
      vector_t v(0);
      if ( x[0] < 0.0 ) {
        d = 1.0;
        p = 1.0;
      }
      else {
        d = 0.125;
        p = 0.1;
      }    
      return std::make_tuple( d, v, p );
    };

  // setup an equation of state
  eos_t eos( /* gamma */ 1.4, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( 
    num_cells_x, num_cells_y, num_cells_z,
    length_x, length_y, length_z );

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // create a single boundary right now
  auto bc_type = boundary_condition_t< mesh_t::num_dimensions >();
  
  // install each boundary
  // -  both +ve and -ve side boundaries can be installed at once since 
  //    they will never overlap

  // x-boundaries
  auto bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[0] == -0.5*length_x || fx[0] == 0.5*length_x );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // y-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[1] == -0.5*length_y || fx[1] == 0.5*length_y );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // z-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[2] == -0.5*length_z || fx[2] == 0.5*length_z );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );

  //===========================================================================
  // SEDOV Inputs
  //===========================================================================

#elif defined(SEDOV_2D)

  using mesh_t = mesh_2d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix  = "sedov_2d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 100;

  // the grid dimensions
  constexpr size_t num_cells_x = 30;
  constexpr size_t num_cells_y = 30;

  constexpr real_t length_x = 1.2;
  constexpr real_t length_y = 1.2;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 0.25, .volume = 0.1, .growth = 1.01 };
  constexpr real_t initial_time_step = 1.e-5;
  constexpr real_t final_time = 1.0;
  size_t max_steps = 1e6;

  // the value of gamma
  constexpr real_t gamma = 1.4;

  // some length scales
  auto dx = length_x / num_cells_x;
  auto dy = length_y / num_cells_y;

  // compute a reference volume
  auto vol = dx * dy;

  // compute a radial size
  auto delta_r = std::numeric_limits<real_t>::epsilon() + 
    std::sqrt( math::sqr( dx/2 ) + math::sqr( dy/2 );
  
  // this is a lambda function to set the initial conditions
  auto ics = [&delta_r,&vol] ( const auto & x )
    {
      constexpr real_t e0 = 0.244816;
      real_t d = 1.0;
      vector_t v = 0;
      real_t p = 1.e-6;
      auto r = sqrt( x[0]*x[0] + x[1]*x[1] );
      if ( r < delta_r  ) 
        p = (gamma - 1) * d * e0 / vol;
      return std::make_tuple( d, v, p );
    };
  
  // setup an equation of state
  eos_t eos(  gamma, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( 
    num_cells_x, num_cells_y, 0.0, 0.0, length_x, length_y 
  );

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // create a single boundary right now
  auto bc_type = boundary_condition_t< mesh_t::num_dimensions >();
  
  // install each boundary
  // -  both +ve and -ve side boundaries can be installed at once since 
  //    they will never overlap

  // x-boundaries
  auto bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[0] == 0.0 || fx[0] == length_x );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // y-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[1] == 0.0 || fx[1] == length_y );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );


  //===========================================================================
  // SEDOV Inputs
  //===========================================================================

#elif defined(SEDOV_3D)

  using mesh_t = mesh_3d_t;
  using real_t = mesh_t::real_t;
  using vector_t = mesh_t::vector_t;

  // the case prefix
  std::string prefix  = "sedov_3d";
  std::string postfix = "vtk";

  // output frequency
  size_t output_freq = 10;

  // the grid dimensions
  constexpr size_t num_cells_x = 20;
  constexpr size_t num_cells_y = 20;
  constexpr size_t num_cells_z = 20;

  constexpr real_t length_x = 1.2;
  constexpr real_t length_y = 1.2;
  constexpr real_t length_z = 1.2;
  
  // the CFL and final solution time
  constexpr time_constants_t 
    CFL = { .accoustic = 0.25, .volume = 0.1, .growth = 1.01 };
  constexpr real_t initial_time_step = 1.e-5;
  constexpr real_t final_time = 1.0;
  size_t max_steps = 1e6;

  // the value of gamma
  constexpr real_t gamma = 1.4;
  
  // some length scales
  auto dx = length_x / num_cells_x;
  auto dy = length_y / num_cells_y;
  auto dz = length_z / num_cells_z;

  // compute a reference volume
  auto vol = dx * dy * dz;

  // compute a radial size
  auto delta_r = std::numeric_limits<real_t>::epsilon() + 
    std::sqrt( math::sqr( dx/2 ) + math::sqr( dy/2 ) + math::sqr( dz/2 ) );
  
  // this is a lambda function to set the initial conditions
  auto ics = [&delta_r,&vol] ( const auto & x )
    {
      constexpr real_t e0 = 0.106384;
      real_t d = 1.0;
      vector_t v = 0;
      real_t p = 1.e-6;
      auto r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
      if ( r < delta_r  ) 
        p = (gamma - 1) * d * e0 / vol;
      return std::make_tuple( d, v, p );
    };
  
  // setup an equation of state
  eos_t eos(  gamma, /* cv */ 1.0 ); 

  // this is the mesh object
  auto mesh = mesh::box<mesh_t>( 
    num_cells_x, num_cells_y, num_cells_z, 0.0, 0.0, 0.0, 
    length_x, length_y, length_z 
  );

  // the boundary mapper
  boundary_map_t< mesh_t::num_dimensions > boundaries;

  // create a single boundary right now
  auto bc_type = boundary_condition_t< mesh_t::num_dimensions >();
  
  // install each boundary
  // -  both +ve and -ve side boundaries can be installed at once since 
  //    they will never overlap

  // x-boundaries
  auto bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[0] == 0.0 || fx[0] == length_x );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // y-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[1] == 0.0 || fx[1] == length_y );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );
  // z-boundaries
  bc_key = mesh.install_boundary( 
    [](auto f) { 
      if ( f->is_boundary() ) {
        const auto & fx = f->midpoint();
        return ( fx[2] == 0.0 || fx[2] == length_z );
      }
      return false; 
    } 
  );
  boundaries.emplace( bc_key, &bc_type );


#else

#  pragma message("NO CASE DEFINED!")

#endif

  //===========================================================================
  // Parse arguments
  //===========================================================================

  // the usage stagement
  auto print_usage = [&argv]() {
    std::cout << "Usage: " << argv[0] 
              << " [--max_steps MAX_STEPS]"
              << " [--output_freq OUTPUT_FREQ]"
              << " [--case CASE_NAME]"
              << " [--postfix EXTENSION]"
              << std::endl << std::endl;
    std::cout << "\t--max_steps MAX_STEPS:\tOverride the number of "
              << "iterations with MAX_STEPS." << std::endl;
    std::cout << "\t--output_freq OUTPUT_FREQ:\tOverride the frequency "
              << "of output to occur ever OUTPUT_FREQ." << std::endl;
    std::cout << "\t--case CASE_NAME:\tOverride the case name "
              << "with CASE_NAME." << std::endl;
    std::cout << "\t--postfix EXTENSION:\tOverride the output extension "
              << "with EXTENSION." << std::endl;
  };

  while (1) {
    
    // Define the options
    static struct option long_options[] =
      {
        {"help",              no_argument, 0, 'h'},
        {"max_steps",   required_argument, 0, 'n'},
        {"output_freq", required_argument, 0, 'f'},
        {"case",        required_argument, 0, 'c'},
        {"postfix",     required_argument, 0, 'p'},
        {0, 0, 0, 0}
      };
      // getopt_long stores the option index here.
      int option_index = 0;

      auto c = getopt_long (argc, argv, "h", long_options, &option_index);

      // Detect the end of the options.
      if (c == -1) break;

      switch (c) {
        case 'c':
          short_warning ("Overriding case name with \"" << optarg << "\"" );
          prefix = optarg;
          break;

        case 'p':
          short_warning ("Overriding extension with \"" << optarg << "\"" );
          postfix = optarg;
          break;

        case 'f':
          short_warning ("Overriding output frequency to \"" << optarg << "\"" );
          output_freq = atoi( optarg );
          break;

        case 'n':
          short_warning ("Overriding max iterations to \"" << optarg << "\"" );
          max_steps = atoi( optarg );
          break;

        case 'h':
          print_usage();
          return 0;

        case '?':
          // getopt_long already printed an error message.
          print_usage();
          return 1;

        default:
          abort ();
        }
    }

  //===========================================================================
  // Mesh Setup
  //===========================================================================  

  // now check the mesh
  mesh.is_valid();
  cout << mesh;
  
  //===========================================================================
  // Field Creation
  //===========================================================================

  // start the timer
  auto tstart = utils::get_wall_time();

  // type aliases
  using matrix_t = matrix_t< mesh_t::num_dimensions >;
  using flux_data_t = flux_data_t< mesh_t::num_dimensions >;

  // create some field data.  Fields are registered as struct of arrays.
  // this allows us to access the data in different patterns.
  register_data(mesh, hydro, cell_mass,       real_t, dense, 1, cells);
  register_data(mesh, hydro, cell_pressure,   real_t, dense, 1, cells);
  register_data(mesh, hydro, cell_velocity, vector_t, dense, 2, cells);

  register_data(mesh, hydro, cell_density,         real_t, dense, 1, cells);
  register_data(mesh, hydro, cell_internal_energy, real_t, dense, 2, cells);
  register_data(mesh, hydro, cell_temperature,     real_t, dense, 1, cells);
  register_data(mesh, hydro, cell_sound_speed,     real_t, dense, 1, cells);

  register_data(mesh, hydro, cell_residual, flux_data_t, dense, 1, cells);

  register_data(mesh, hydro, node_coordinates, vector_t, dense, 1, vertices);
  register_data(mesh, hydro, node_velocity, vector_t, dense, 1, vertices);
  
  register_data(mesh, hydro, corner_matrix, matrix_t, dense, 1, corners);
  register_data(mesh, hydro, corner_normal, vector_t, dense, 1, corners);

  // register the time step and set a cfl
  register_data( mesh, hydro, time_step, real_t, global, 1 );
  register_data( mesh, hydro, cfl, time_constants_t, global, 1 );
  
  *get_accessor( mesh, hydro, time_step, real_t, global, 0 ) = initial_time_step;
  *get_accessor( mesh, hydro, cfl, time_constants_t, global, 0 ) = CFL;


  // register state a global eos
  register_data( mesh, hydro, eos, eos_t, global, 1 );
  *get_accessor( mesh, hydro, eos, eos_t, global, 0 ) = eos;

  // set the persistent variables, i.e. the ones that will be plotted
  get_accessor(mesh, hydro, cell_mass,       real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, cell_pressure,   real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, cell_velocity, vector_t, dense, 0).attributes().set(persistent);

  get_accessor(mesh, hydro, cell_density,         real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, cell_internal_energy, real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, cell_temperature,     real_t, dense, 0).attributes().set(persistent);
  get_accessor(mesh, hydro, cell_sound_speed,     real_t, dense, 0).attributes().set(persistent);

  get_accessor(mesh, hydro, node_velocity, vector_t, dense, 0).attributes().set(persistent);
  

  //===========================================================================
  // Initial conditions
  //===========================================================================
  
  // now call the main task to set the ics.  Here we set primitive/physical 
  // quanties
  apps::hydro::initial_conditions( mesh, ics );
  

  // Update the EOS
  apps::hydro::update_state_from_pressure( mesh );

  //===========================================================================
  // Pre-processing
  //===========================================================================



  // now output the solution
  if ( output_freq > 0 )
    apps::hydro::output(mesh, prefix, postfix, 1);
  
  //===========================================================================
  // Residual Evaluation
  //===========================================================================


  // get the current time
  auto soln_time = mesh.time();
  auto time_cnt  = mesh.time_step_counter();

  // get the time step accessor
  auto time_step = get_accessor( mesh, hydro, time_step, real_t, global, 0 );   

  // a counter for this session
  size_t num_steps = 0; 

  for ( ; (num_steps < max_steps && soln_time < final_time); ++num_steps ) 
  {   

    //--------------------------------------------------------------------------
    // Begin Time step
    //--------------------------------------------------------------------------

    // Save solution at n=0
    apps::hydro::save_coordinates( mesh );
    apps::hydro::save_solution( mesh );

    //--------------------------------------------------------------------------
    // Predictor step : Evaluate Forces at n=0
    //--------------------------------------------------------------------------

    // estimate the nodal velocity at n=0
    apps::hydro::estimate_nodal_state( mesh );

    // evaluate corner matrices and normals at n=0 
    apps::hydro::evaluate_corner_coef( mesh );

    // compute the nodal velocity at n=0
    apps::hydro::evaluate_nodal_state( mesh, boundaries );

    // compute the fluxes
    apps::hydro::evaluate_forces( mesh );

    //--------------------------------------------------------------------------
    // Time step evaluation
    //--------------------------------------------------------------------------

    // compute the time step
    std::string limit_string;
    apps::hydro::evaluate_time_step( mesh, limit_string );
    
    // access the computed time step and make sure its not too large
    *time_step = std::min( *time_step, final_time - soln_time );       

    cout << std::string(60, '=') << endl;
    auto ss = cout.precision();
    cout.setf( std::ios::scientific );
    cout.precision(6);
    cout << "| " << std::setw(8) << "Step:"
         << " | " << std::setw(13) << "Time:"
         << " | " << std::setw(13) << "Step Size:"
         << " | " << std::setw(13) << "Limit:"
         << " |" << std::endl;
    cout << "| " << std::setw(8) << time_cnt+1
         << " | " << std::setw(13) << soln_time + (*time_step)
         << " | " << std::setw(13) << *time_step
         << " | " << std::setw(13) << limit_string
         << " |" << std::endl;
    cout.unsetf( std::ios::scientific );
    cout.precision(ss);

// #define USE_FIRST_ORDER_TIME_STEPPING
#ifndef USE_FIRST_ORDER_TIME_STEPPING // set to 0 for first order

    //--------------------------------------------------------------------------
    // Move to n+1/2
    //--------------------------------------------------------------------------

    // move the mesh to n+1/2
    apps::hydro::move_mesh( mesh, 0.5 );

    // update solution to n+1/2
    apps::hydro::apply_update( mesh, 0.5, true );

    // Update derived solution quantities
    apps::hydro::update_state_from_energy( mesh );

    //--------------------------------------------------------------------------
    // Corrector : Evaluate Forces at n=1/2
    //--------------------------------------------------------------------------

    // evaluate corner matrices and normals at n=1/2
    apps::hydro::evaluate_corner_coef( mesh );

    // compute the nodal velocity at n=1/2
    apps::hydro::evaluate_nodal_state( mesh, boundaries );

    // compute the fluxes
    apps::hydro::evaluate_forces( mesh );

    //--------------------------------------------------------------------------
    // Move to n+1
    //--------------------------------------------------------------------------

    // restore the solution to n=0
    apps::hydro::restore_coordinates( mesh );
    apps::hydro::restore_solution( mesh );

#endif

    // move the mesh to n+1
    apps::hydro::move_mesh( mesh, 1.0 );
    
    // update solution to n+1
    apps::hydro::apply_update( mesh, 1.0, false );

    // Update derived solution quantities
    apps::hydro::update_state_from_energy( mesh );

    //--------------------------------------------------------------------------
    // End Time step
    //--------------------------------------------------------------------------

    // update time
    soln_time = mesh.increment_time( *time_step );
    time_cnt = mesh.increment_time_step_counter();
  
    // now output the solution
    apps::hydro::output(mesh, prefix, postfix, output_freq);

  }


  //===========================================================================
  // Post-process
  //===========================================================================
    
  // now output the solution
  if ( (output_freq > 0) && (time_cnt%output_freq != 0) )
    apps::hydro::output(mesh, prefix, postfix, 1);

  cout << "Final solution time is " 
       << std::scientific << std::setprecision(6) << soln_time
       << " after " << num_steps << " steps." << std::endl;

  
  auto tdelta = utils::get_wall_time() - tstart;
  std::cout << "Elapsed wall time is " << std::setprecision(4) << std::fixed 
            << tdelta << "s." << std::endl;
  
  // success
  return 0;

}

