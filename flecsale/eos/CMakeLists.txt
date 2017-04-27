#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(eos_HEADERS
  eos_base.h
  ideal_gas.h
)


cinch_install_headers(
  FILES ${eos_HEADERS}
  DESTINATION include/flecsale/eos)


cinch_add_unit( test_eos
  SOURCES test/ideal_gas.cc
)
