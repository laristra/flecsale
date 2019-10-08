/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Utilities to process command line arguments.
///////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include <flecsi-sp/utils/parse_arguments.h>

// system includes
#include <getopt.h>

namespace apps {
namespace common {

///////////////////////////////////////////////////////////////////////////////
//! \brief Process the argument list for this app.
///////////////////////////////////////////////////////////////////////////////
inline auto process_arguments(int argc, char** argv)
{

  // the usage stagement
  auto print_usage = [&argv]() {
    std::cout << "Usage: " << argv[0]
              << " [--input INPUT_FILE]"
              << " [--max_entries MAX_ENTRIES]"
              << " [--mesh MESH_FILE]"
              << " [--help]"
              << std::endl << std::endl;
    std::cout << "\t--input_file INPUT_FILE:\t Override the input file "
              << "with INPUT_FILE." << std::endl;
    std::cout << "\t--max_entries MAX_ENTRIES:\t Override the maximum number "
              << "of sparse entries per entity with ENTRIES." << std::endl;
    std::cout << "\t--mesh MESH_FILE:\t Override the mesh file "
              << "with MESH_FILE." << std::endl;
    std::cout << "\t--help:\t Print a help message." << std::endl;
  };

  // Define the options
  struct option long_options[] =
    {
      {"help",            no_argument, 0, 'h'},
      {"input_file",    required_argument, 0, 'f'},
      {"max_entries",   required_argument, 0, 'e'},
      {"mesh",      required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };
  const char * short_options = "hf:e:m:";

  // parse the arguments
  auto args =
    flecsi_sp::utils::parse_arguments(argc, argv, long_options, short_options);

  // process the simple ones
  if ( args.count("h") ) print_usage();

  return args;
}

} // namespace
} // namespace
