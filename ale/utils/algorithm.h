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
 * \file type_traits.h
 * 
 * \brief Extensions to C++11 that are part of C++14.
 *
 ******************************************************************************/
#pragma once

//! system includes
#include <alorithm>

namespace ale {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! Stuff not in C++
////////////////////////////////////////////////////////////////////////////////

//==============================================================================
// perform a set intersection, but this one sorts first
//==============================================================================
template<class InputIt1, class InputIt2, class OutputIt>
OutputIt unsorted_set_intersection( InputIt1 first1, InputIt1 last1,
                                    InputIt2 first2, InputIt2 last2,
                                    OutputIt d_first )
{
  using value_type1 = std::iterator_traits<InputIt1>::value_type;
  using value_type2 = std::iterator_traits<InputIt2>::value_type;
  
  std::vector< value_type1 > sorted1( first1, last1 );
  std::vector< value_type2 > sorted2( first2, last2 );
  
  std::sort( sorted1.begin(), sorted1.end() );
  std::sort( sorted2.begin(), sorted2.end() );

  return std::set_intersection( sorted1.begin(), sorted1.end(), 
                                sorted2.begin(), sorted2.end(),
                                d_first );
}

} // namespace
} // namespace
