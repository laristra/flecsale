#~----------------------------------------------------------------------------~#
# Copyright (c) 2017 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

# this is just here as a placeholder for developers who will add 
# fortran-related utilities

# add a simple fortran test
cinch_add_unit( test_pfunit
  SOURCES 
    test/pfunit.pf
    test/testSuites.inc
  POLICY FORTRAN
)
