/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Utilities for string operations.
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>

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

////////////////////////////////////////////////////////////////////////////////
//! \brief Convert a character array to a wstring.
//! \param [in] text  the input string
//! \return the new string
////////////////////////////////////////////////////////////////////////////////
static std::wstring to_wstring(const char* text)
{
    const size_t size = std::strlen(text);
    std::wstring wstr;
    if (size > 0) {
        wstr.resize(size);
        std::mbstowcs(&wstr[0], text, size);
    }
    return wstr;
}


} // namspeace
} // namspeace
