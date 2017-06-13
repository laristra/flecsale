#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro_3d
  $<TARGET_OBJECTS:common> 
  inputs.cc
  driver.cc
  ${FleCSI_RUNTIME_DRIVER}
  ${FleCSI_RUNTIME_MAIN}
)
target_link_libraries( hydro_3d flecsale )

if (LUA_FOUND)

  create_regression_test( 
    NAME shock_box_3d
    COMMAND $<TARGET_FILE:hydro_3d> -f ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d.lua
    COMPARE shock_box_3d0000011.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000011.dat.std 
  )
  
  create_regression_test( 
    NAME shock_box_3d_omp4
    COMMAND $<TARGET_FILE:hydro_3d> -f ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d.lua
    THREADS 4
    COMPARE shock_box_3d0000011.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000011.dat.std 
  )

else()

  create_regression_test( 
    NAME shock_box_3d
    COMMAND $<TARGET_FILE:hydro_3d> 
    COMPARE shock_box_3d0000011.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000011.dat.std 
  )
  
  create_regression_test( 
    NAME shock_box_3d_omp4
    COMMAND $<TARGET_FILE:hydro_3d> 
    THREADS 4
    COMPARE shock_box_3d0000011.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d0000011.dat.std 
  )

endif()
