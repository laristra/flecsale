#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( maire_hydro_2d
  inputs.cc
  driver.cc
  ${FLECSI_RUNTIME_DRIVER}
  ${FLECSI_RUNTIME_MAIN}
)
target_link_libraries( maire_hydro_2d common ${ALE_LIBRARIES} )

if (LUA_FOUND) 

  create_regression_test( 
    NAME sedov_maire_2d
    COMMAND $<TARGET_FILE:maire_hydro_2d> -f ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d.lua
    COMPARE sedov_2d0000020.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d0000020.dat.std 
  )
  
  create_regression_test( 
    NAME sedov_maire_2d_omp4
    COMMAND $<TARGET_FILE:maire_hydro_2d> -f ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d.lua
    THREADS 4
    COMPARE sedov_2d0000020.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d0000020.dat.std 
  )

else()

  create_regression_test( 
    NAME sedov_maire_2d
    COMMAND $<TARGET_FILE:maire_hydro_2d>
    COMPARE sedov_2d0000020.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d0000020.dat.std 
  )
  
  create_regression_test( 
    NAME sedov_maire_2d_omp4
    COMMAND $<TARGET_FILE:maire_hydro_2d>
    THREADS 4
    COMPARE sedov_2d0000020.dat 
    STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d0000020.dat.std 
  )

endif()
