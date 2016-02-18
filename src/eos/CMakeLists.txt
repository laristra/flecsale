#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(eos_HEADERS
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

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
