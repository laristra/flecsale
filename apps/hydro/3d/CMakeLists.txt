#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set_property(SOURCE driver.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE inputs.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../specializations/flecsi-sp/burton/burton_specialization_init.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../flecsi/flecsi/execution/legion/runtime_driver.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../flecsi/flecsi/execution/legion/runtime_main.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
add_executable( hydro_3d
  $<TARGET_OBJECTS:apps_common> 
  driver.cc
  inputs.cc
  ${FLECSI_SP_BURTON_SPECIALIZATION_INIT}
  ${FLECSALE_RUNTIME_DRIVER}
  ${FLECSALE_RUNTIME_MAIN}
)
target_link_libraries( hydro_3d FleCSALE -lkokkoscore )
target_compile_definitions( 
  hydro_3d
  PRIVATE
    FLECSI_SP_BURTON_MESH_DIMENSION=3
    FLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
)

if ( FLECSALE_ENABLE_LUA )
  create_regression_test(
    NAME flecsale_shock_box_3d_1procs_with_input
    COMMAND
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:hydro_3d> ${MPIEXEC_POSTFLAGS} -m
      ${FLECSALE_DATA_DIR}/meshes/sedov_20x20x20.g
      -f ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d.lua
    COMPARE
      shock_box_3d-solution_rank000000.txt
    STANDARD
      ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d-solution_rank000000.txt.1procs.std
  )
endif()

create_regression_test(
  NAME flecsale_shock_box_3d_1procs
  COMMAND
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
  $<TARGET_FILE:hydro_3d> ${MPIEXEC_POSTFLAGS} -m
  ${FLECSALE_DATA_DIR}/meshes/sedov_20x20x20.g
  COMPARE
  shock_box_3d-solution_rank000000.txt
  STANDARD
  ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d-solution_rank000000.txt.1procs.std
  )

create_regression_test(
  NAME flecsale_shock_box_3d_2procs
  COMMAND
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
  $<TARGET_FILE:hydro_3d> ${MPIEXEC_POSTFLAGS} -m
  ${FLECSALE_DATA_DIR}/meshes/sedov_20x20x20.g
  COMPARE
  shock_box_3d-solution_rank000000.txt
  shock_box_3d-solution_rank000001.txt
  STANDARD
  ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d-solution_rank000000.txt.2procs.std
  ${CMAKE_CURRENT_SOURCE_DIR}/shock_box_3d-solution_rank000001.txt.2procs.std
  )
