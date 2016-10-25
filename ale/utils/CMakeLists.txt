#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(utils_HEADERS
  algorithm.h
  array_ref.h
  array_view.h
  const_string.h
  exceptions.h
  errors.h
  filter_iterator.h
  fixed_vector.h
  string_utils.h
  static_for.h
  tasks.h
  template_helpers.h detail/template_helpers_impl.h
  time_utils.h
  tuple_for_each.h   detail/tuple_for_each_impl.h
  tuple_visit.h      detail/tuple_visit_impl.h
  tuple_zip.h        detail/tuple_zip_impl.h
  type_traits.h      detail/type_traits_impl.h
  zip.h              detail/zip_impl.h
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
