#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(eos_HEADERS
  eos_base.h
  ideal_gas.h
  
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)


cinch_add_unit( test_eos
  SOURCES test/ideal_gas.cc
)
