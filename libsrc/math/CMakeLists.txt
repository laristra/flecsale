#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(math_HEADERS
  array.h
  constants.h
  general.h  detail/general_impl.h
  matrix.h
  multi_array.h
  tuple.h
  vector.h
  PARENT_SCOPE
)

cinch_add_unit(math
    SOURCES 
      test/tuple.cc 
      test/vector.cc
      test/matrix.cc
)
