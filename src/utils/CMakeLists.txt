#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(utils_HEADERS
  check_types.h
  static_for_each.h
  tasks.h
  template_helpers.h
  tuple.h
  tuple_zip.h
  vector.h
  zip.h
  PARENT_SCOPE
)

set(utils_SOURCES
  dummy.c
  PARENT_SCOPE
)

cinch_add_unit(vector
    SOURCES test/vector.cc
)

cinch_add_unit(tuple
    SOURCES test/tuple.cc
)

cinch_add_unit(zip_iterator
    SOURCES test/zip.cc
)


cinch_add_unit(tasks
    SOURCES test/tasks.cc
)


#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
