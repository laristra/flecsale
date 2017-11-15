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
  functional.h
  lua_utils.h
  python_utils.h
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
  
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

cinch_add_unit(test_utils
    SOURCES 
      test/array_view.cc
      test/caliper.cc
      test/fixed_vector.cc
      test/lua_utils.cc
      test/python_utils.cc
      test/static_for.cc
      test/tasks.cc
      test/tuple_for_each.cc
      test/tuple_visit.cc
      test/tuple_zip.cc
      test/zip.cc
    INPUTS 
      test/lua_test.lua
      test/python_test.py
)

# make sure the test was actually made
if (TARGET test_utils)
  
  # now link the necessary libraries
  if (PYTHON_LIBRARIES) 
    target_link_libraries(test_utils ${PYTHON_LIBRARIES})
  endif()
  if (LUA_LIBRARIES) 
    target_link_libraries(test_utils ${LUA_LIBRARIES})
  endif()
  if (Caliper_LIBRARIES) 
    target_link_libraries(test_utils ${Caliper_LIBRARIES})
  endif()

endif()

