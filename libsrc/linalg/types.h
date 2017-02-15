/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Defines the types used for the linear algebra solvers.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "utils/array_view.h"


namespace ale {
namespace linalg {

//! \brief The matrix type for all linear solvers in this namespace.
template< typename T >
using matrix_view = utils::array_view<T,2>;

//! \brief The vector type for all linear solvers in this namespace.
template< typename T >
using vector_view = utils::array_view<T,1>;

} // namespace
} // namespace


