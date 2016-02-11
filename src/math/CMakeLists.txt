#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(math_HEADERS
  array.h
  operators.h  detail/operators.h
  tuple.h
  vector.h
  PARENT_SCOPE
)

cinch_add_unit(test_math
    SOURCES 
      test/tuple.cc 
      test/vector.cc
      test/matrix.cc
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
