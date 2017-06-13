#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( maire_hydro_3d
  $<TARGET_OBJECTS:common> 
  inputs.cc
  driver.cc
  ${FleCSI_RUNTIME_DRIVER}
  ${FleCSI_RUNTIME_MAIN}
)
target_link_libraries( maire_hydro_3d flecsale )

if (LUA_FOUND) 

  create_regression_test( 
    NAME sedov_maire_3d
    COMMAND $<TARGET_FILE:maire_hydro_3d> -f ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d.lua
    COMPARE sedov_3d0000010.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
  )
  
  create_regression_test( 
    NAME sedov_maire_3d_omp4
    COMMAND $<TARGET_FILE:maire_hydro_3d> -f ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d.lua
    THREADS 4
    COMPARE sedov_3d0000010.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
  )

else()

  create_regression_test( 
    NAME sedov_maire_3d
    COMMAND $<TARGET_FILE:maire_hydro_3d> 
    COMPARE sedov_3d0000010.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
  )
  
  create_regression_test( 
    NAME sedov_maire_3d_omp4
    COMMAND $<TARGET_FILE:maire_hydro_3d> 
    THREADS 4
    COMPARE sedov_3d0000010.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_3d0000010.dat.std 
  )

endif()
