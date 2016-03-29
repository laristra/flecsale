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
 * \file string_utils.h
 * 
 * \brief Utilities for string operations.
 *
 ******************************************************************************/
#pragma once

namespace ale {
namespace utils {


////////////////////////////////////////////////////////////////////////////////
//! \brief replace all occuraces of "from" to "to"
//! \param [in] str  the input string
//! \param [in] from the string to search for
//! \param [in] to   the string to replace "from" with
//! \return the new string
////////////////////////////////////////////////////////////////////////////////
inline
auto replace_all(std::string str, const std::string & from, const std::string & to) {
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}

} // namspeace
} // namspeace
