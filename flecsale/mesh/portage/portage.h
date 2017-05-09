/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief The main include for using portage with the meshes defined in this
///   library.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsale/mesh/portage/portage_mesh.h"
#include "flecsale/mesh/portage/portage_state.h"

// portage library includes
#define PORTAGE_SERIAL_ONLY
#include <portage/driver/driver.h>
#undef PORTAGE_SERIAL_ONLY
