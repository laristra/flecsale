#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(eqns_HEADERS
  euler_eqns.h
  flux.h
  lagrange_eqns.h
  
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

cinch_add_unit( test_eqns
  SOURCES 
    test/euler_eqns.cc
)
