#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(math_HEADERS
  tuple.h
  vector.h
  PARENT_SCOPE
)

cinch_add_unit(vector
    SOURCES test/vector.cc
)

cinch_add_unit(tuple
    SOURCES test/tuple.cc
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
