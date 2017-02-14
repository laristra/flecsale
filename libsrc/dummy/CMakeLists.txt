#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

# these are all the headers within this sub-package/sub-folder
set(dummy_HEADERS
  dummy.h
  PARENT_SCOPE
)


# these are all the sources within this sub-package/sub-folder
set(dummy_SOURCES
  dummy.cc
  PARENT_SCOPE
)

# if the test uses a library, uncomment the line below and replace
# "some_library" with the required library

# set( dummy_LIBS some_library )

# add a unit test.  test code for dummy.h/dummy.cc goes in a folder
# called test inside this directory
cinch_add_unit(dummy
    SOURCES 
      test/dummy.cc
    LIBRARIES ${ALE_LIBRARIES}
)

