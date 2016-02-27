#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(math_HEADERS
  array/array.h
  array/multi_array.h
  array/super_class.h
  math.h    detail/math.h
  tuple.h
  vector.h
  matrix.h
  PARENT_SCOPE
)

cinch_add_unit(math
    SOURCES 
      test/tuple.cc 
      test/vector.cc
      test/matrix.cc
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
