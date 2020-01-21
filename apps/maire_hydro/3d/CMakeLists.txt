#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set_property(SOURCE driver.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE inputs.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../specializations/flecsi-sp/burton/burton_specialization_init.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../flecsi/flecsi/execution/legion/runtime_driver.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../flecsi/flecsi/execution/legion/runtime_main.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
add_executable( maire_hydro_3d
  $<TARGET_OBJECTS:apps_common> 
  driver.cc
  inputs.cc
  ${FLECSI_SP_BURTON_SPECIALIZATION_INIT}
  ${FLECSALE_RUNTIME_DRIVER}
  ${FLECSALE_RUNTIME_MAIN}
)
target_link_libraries( maire_hydro_3d FleCSALE )
target_compile_definitions( 
  maire_hydro_3d
  PRIVATE
    FLECSI_SP_BURTON_MESH_DIMENSION=3
    FLECSI_SP_BURTON_MESH_EXTRAS
    FLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
)

add_test( 
  NAME flecsale_sedov_3d_1procs
  COMMAND mpirun -n 1 $<TARGET_FILE:maire_hydro_3d> -m ${FLECSALE_DATA_DIR}/meshes/sedov_20x20x20.g
  #COMPARE shock_box_2d0000007.dat 
  #STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_2d0000007.dat.std 
)

add_test( 
  NAME flecsale_sedov_3d_2procs
  COMMAND mpirun -n 2 $<TARGET_FILE:maire_hydro_3d> -m ${FLECSALE_DATA_DIR}/meshes/sedov_20x20x20.g
  #COMPARE shock_box_2d0000007.dat 
  #STANDARD ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_2d0000007.dat.std 
)
