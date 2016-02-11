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
 * \file operators.h
 * 
 * \brief Some helper functions for template foo magic.
 *
 ******************************************************************************/
#pragma once

namespace ale {
namespace math {
  namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Template helper to compute averages
////////////////////////////////////////////////////////////////////////////////

//! \brief average operator.
template< class T >
constexpr void average( T & res )
{ 
  // nothing left to do
}

//! \brief average operator.
template< class T, class U, class ... Args >
constexpr void average( T & res, U && u, Args&&... args )
{ 
  res += u;
  average(res, std::forward<Args>(args)...); 
}


}  // namespace
} // namespace
} // namespace



/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
