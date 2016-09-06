#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(eos_HEADERS
  euler_eqns.h
  flux.h
  lagrange_eqns.h
  PARENT_SCOPE
)

#set(eos_SOURCES
#  dummy.c
#  PARENT_SCOPE
#)

cinch_add_unit( eqns
  SOURCES 
    test/euler_eqns.cc
)
