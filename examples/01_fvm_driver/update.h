/*~--------------------------------------------------------------------------~*
 *     _   ______________     ___    __    ______
 *    / | / / ____/ ____/    /   |  / /   / ____/
 *   /  |/ / / __/ /  ______/ /| | / /   / __/   
 *  / /|  / /_/ / /__/_____/ ___ |/ /___/ /___   
 * /_/ |_/\____/\____/    /_/  |_/_____/_____/   
 * 
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef update_h
#define update_h

#include "types.h"

int32_t update(burton_mesh_t & mesh, const accessor_t & p,
  const accessor_t & d, accessor_t & r) {

  for(auto z: mesh.cells()) {
    r[z] = p[z] + d[z];
  } // for

  return 0;
} // driver

#endif // update_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
