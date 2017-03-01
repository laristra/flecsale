#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(eos_HEADERS
  eos_base.h
  ideal_gas.h
  PARENT_SCOPE
)

#set(eos_SOURCES
#  dummy.c
#  PARENT_SCOPE
#)

cinch_add_unit( eos
  SOURCES test/ideal_gas.cc
)
