#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(utils_HEADERS
  check_types.h    detail/check_types.h
  const_string.h
  exceptions.h
  errors.h
  static_for.h     detail/static_for.h
  tasks.h          detail/tasks.h
  template_helpers.h  detail/template_helpers.h
  tuple_for_each.h detail/tuple_for_each.h
  tuple_visit.h    detail/tuple_visit.h
  tuple_zip.h      detail/tuple_zip.h
  write_binary.h
  zip.h            detail/zip.h
  PARENT_SCOPE
)

set(utils_SOURCES
  vtk.cc
  PARENT_SCOPE
)

cinch_add_unit(utils
    SOURCES 
      test/static_for.cc
      test/tasks.cc
      test/tuple_for_each.cc
      test/tuple_visit.cc
      test/tuple_zip.cc
      test/zip.cc
)


#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
