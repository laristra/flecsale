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

namespace apps {
namespace common {


///////////////////////////////////////////////////////////////////////////////
//! \brief Parse the argument list
///////////////////////////////////////////////////////////////////////////////
static auto parse_arguments(
  int argc, 
  char** argv, 
  const option * long_options, 
  const char * short_options
) {

  std::map<std::string, std::string> key_value_pair;

  // reset getopts global variable
  optind = 0;

    
    // getopt_long stores the option index here.
  int option_index = 0;

  int c=0;
  int done = 0;

  while (!done && (c != -1)) {
    auto c = 
      getopt_long(argc, argv, short_options, long_options, &option_index);
    auto c_char = static_cast<char>(c);
    auto c_str = flecsale::utils::to_string( c_char );
    switch (c) {
      case (-1):
        break;
      case 0:
        key_value_pair[long_options[option_index].name] =
          optarg ? optarg : "";
      default:
         done = 1;
        break;
    }
    if (optarg) key_value_pair[c_str] = optarg;
    else        key_value_pair[c_str] = "";

    if (optind > argc)
    raise_runtime_error( "Expected argument after options" );

  }

  return key_value_pair;
}

} // namespace
} // namespace
