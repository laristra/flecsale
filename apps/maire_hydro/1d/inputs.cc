/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This is where the compile time inputs are defined.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include "../inputs.h"

#include <flecsi/flecsi/execution/execution.h>
#include <flecsale/eos/ideal_gas.h>

// system includes
#include <limits>

namespace apps {
namespace hydro {

//=============================================================================
// type aliases confined to this translation unit
//=============================================================================
using size_t = inputs_t::size_t;
using real_t = inputs_t::real_t;
using vector_t = inputs_t::vector_t;
using string = std::string;
using eos_t = inputs_t::eos_t;

//=============================================================================
// These constants are not part of the input class, they are globally scoped
// within this translation unit only.
//=============================================================================

// the value of gamma
constexpr real_t gamma = 1.4;

// utility to see if value near another, within tolerance
auto near( real_t a, real_t b, real_t tol = flecsale::config::test_tolerance  )
{ return std::abs(a - b) <= tol; }

//=============================================================================
// Now set the inputs.
//=============================================================================

// the case prefix
string inputs_t::prefix = "sedov_1d";
string inputs_t::postfix = "dat";

// output frequency
size_t inputs_t::output_freq = 20;

// the CFL and final solution time
time_constants_t inputs_t::CFL = 
{ .accoustic = 0.25, .volume = 0.1, .growth = 1.01 };

real_t inputs_t::final_time = 1.0;
real_t inputs_t::initial_time_step = 1.e-5;
size_t inputs_t::max_steps = 20;

// the equation of state
eos_t inputs_t::eos = 
  flecsale::eos::ideal_gas_t<real_t>( 
    /* gamma */ 1.4, /* cv */ 1.0 
  ); 

// this is a static function to set the initial conditions
inputs_t::ics_return_t inputs_t::initial_conditions(
	const mesh_t & mesh, size_t local_id, const real_t &)
{
	const auto & x = mesh.cells()[local_id]->centroid();
	const auto & vol = mesh.cells()[local_id]->volume();

	// compute a radial size
	auto delta_r = std::numerical_limits<real_t>::epsilon() +
		std::abs(vol/2);

	constexpr real_t e0 = 0.244816;
	real_t d = 1.0;
	vector_t v = 0;
	real_t p = 1.e-6;
	auto r = abs( x[0] );
	if ( r < delta_r  )
		p = (gamma - 1) * d * e0 / vol;
	return std::make_tuple( d, v, p );
}


// register a global data object for each boundary
static constexpr x_boundary = 0;
flecsi_register_global_object(
	x_boundary,
	boundaries, // arbitrary name
	boundary_condition_t);

void inputs_t::initialize_boundaries() {
	flecsi_initialize_global_object(
		x_boundary,
		boundaries,
		symmetry_boundary_condition_t);
}

// install each boundary
//
// - both +ve and -ve side boundaries can be installed at once since
//   they will never overlap
// - if they did overlap, then you need to define them seperately or else
//   its hard to count the number of different conditions on points or edges.
void inputs_t::boundary_conditions(
	mesh_t & mesh, const real_t &) {
	mesh.install_boundary(
		[=](auto f) {
			if (f->is_boundary()) {
				const auto & n = f->normal();
				return (near(std::abs(n[0]), 1.0));
			}
			else
				return false;
		},
		x_boundary);
}

} // namespace
} // namespace
