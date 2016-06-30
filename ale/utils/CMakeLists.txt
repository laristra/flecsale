#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(utils_HEADERS
  algorithm.h
  array_ref.h
  array_view.h
  type_traits.h    detail/type_traits.h
  const_string.h
  exceptions.h
  errors.h
  fixed_vector.h
  string_utils.h
  static_for.h     detail/static_for.h
  tasks.h
  template_helpers.h  detail/template_helpers.h
  tuple_for_each.h detail/tuple_for_each.h
  tuple_visit.h    detail/tuple_visit.h
  tuple_zip.h      detail/tuple_zip.h
  zip.h            detail/zip.h
  PARENT_SCOPE
)

#set(utils_SOURCES
#  vtk.cc
#  PARENT_SCOPE
#)

cinch_add_unit(utils
    SOURCES 
      test/array_view.cc
      test/fixed_vector.cc
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
