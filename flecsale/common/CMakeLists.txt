#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(common_HEADERS
  types.h
)

cinch_install_headers(
  FILES ${common_HEADERS}
  DESTINATION include/flecsale/common
)

