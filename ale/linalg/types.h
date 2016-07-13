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
 * \file types.h
 * 
 * \brief Defines the types used for the linear algebra solvers.
 *
 ******************************************************************************/
#pragma once

//! user includes
#include "ale/utils/array_view.h"


namespace ale {
namespace linalg {

//! \brief the matrix type for all linear solvers in this namespace
template< typename T >
using matrix_view = utils::array_view<T,2>;

//! \brief the vector type for all linear solvers in this namespace
template< typename T >
using vector_view = utils::array_view<T,1>;

} // namespace
} // namespace


