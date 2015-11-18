#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#set(hydro_HEADERS
#  euler_eqns.h
#  PARENT_SCOPE
#)

#set(hydro_SOURCES
#  dummy.c
#  PARENT_SCOPE
#)


cinch_add_unit(
  hydro
  SOURCES test/hydro.cc
  LIBRARIES ${IO_LIBRARIES}
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
