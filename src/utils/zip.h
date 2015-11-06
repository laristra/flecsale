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
 * \file zip.h
 * 
 * \brief Provide a zip-like iterator for range-based fors.
 *
 ******************************************************************************/
#pragma once


// system includes
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

namespace ale {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! \brief Combine iterators together using boost
////////////////////////////////////////////////////////////////////////////////
template <typename... T>
auto zip(T&&... containers)
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}

} // namespace
} // namespace
