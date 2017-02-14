/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utilities for parsing arguments.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include <flecsale/utils/errors.h>
#include <flecsale/utils/string_utils.h>

// system libraries
#include <cstring>
#include <getopt.h>
#include <map>
#include <string>

///////////////////////////////////////////////////////////////////////////////
//! \brief Parse the argument list
///////////////////////////////////////////////////////////////////////////////
auto parse_arguments(
  int argc, 
  char** argv, 
  const option * long_options, 
  const char * short_options
) {

  std::map<std::string, std::string> key_value_pair;

  while (1) {
    
    // getopt_long stores the option index here.
    int option_index = 0;

    auto c = 
      getopt_long(argc, argv, short_options, long_options, &option_index);
    auto c_char = static_cast<char>(c);
    auto c_str = ale::utils::to_string( c_char );

    // Detect the end of the options.
    if (c == -1) break;

    // no equivalent short option
    if (c == 0) {
      key_value_pair[long_options[option_index].name] = 
        optarg ? optarg : "";
    }
    // if the flag is not found, getopt prints a warning

    // store the value
    if (optarg) key_value_pair[c_str] = optarg;
    else        key_value_pair[c_str] = "";


  }

  return key_value_pair;
}

