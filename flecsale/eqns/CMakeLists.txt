#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(eqns_HEADERS
  euler_eqns.h
  flux.h
  lagrange_eqns.h
)

mcinch_install_headers(
  FILES ${eqns_HEADERS}
  DESTINATION include/flecsale/eqns
)

mcinch_add_unit( test_eqns
  SOURCES 
    test/euler_eqns.cc
)
