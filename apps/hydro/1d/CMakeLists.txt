#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( hydro_1d
  $<TARGET_OBJECTS:apps_common> 
  driver.cc
  inputs.cc
  ${FLECSI_SP_BURTON_SPECIALIZATION_INIT}
  ${FLECSALE_RUNTIME_DRIVER}
  ${FLECSALE_RUNTIME_MAIN}
)
target_link_libraries( hydro_1d FleCSALE )
target_compile_definitions( 
  hydro_1d
  PRIVATE
    FLECSI_SP_BURTON_MESH_DIMENSION=1
    FLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
)

add_test( 
  NAME flecsale_shock_box_1d_1procs
  COMMAND mpirun -n 1 $<TARGET_FILE:hydro_1d> -m not_used
  #COMPARE shock_box_1d0000007.dat 
  #STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_1d0000007.dat.std 
)

add_test( 
  NAME flecsale_shock_box_1d_2procs
  COMMAND mpirun -n 2 $<TARGET_FILE:hydro_1d> -m not_used
  #COMPARE shock_box_1d0000007.dat 
  #STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_1d0000007.dat.std 
)
