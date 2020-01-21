#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set_property(SOURCE driver.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE inputs.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE inputs_regression.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../specializations/flecsi-sp/burton/burton_specialization_init.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../flecsi/flecsi/execution/legion/runtime_driver.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
set_property(SOURCE ../../../flecsi/flecsi/execution/legion/runtime_main.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
add_executable( maire_hydro_2d
  $<TARGET_OBJECTS:apps_common> 
  driver.cc
  inputs.cc
  ${FLECSI_SP_BURTON_SPECIALIZATION_INIT}
  ${FLECSALE_RUNTIME_DRIVER}
  ${FLECSALE_RUNTIME_MAIN}
)
target_link_libraries( maire_hydro_2d FleCSALE )
target_compile_definitions( 
  maire_hydro_2d
  PRIVATE
    FLECSI_SP_BURTON_MESH_DIMENSION=2
    FLECSI_SP_BURTON_MESH_EXTRAS
    FLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
)


add_executable( maire_hydro_2d_regression
  $<TARGET_OBJECTS:apps_common>
  driver.cc
  inputs_regression.cc
  ${FLECSI_SP_BURTON_SPECIALIZATION_INIT}
  ${FLECSALE_RUNTIME_DRIVER}
  ${FLECSALE_RUNTIME_MAIN}
)
target_link_libraries( maire_hydro_2d_regression FleCSALE )
target_compile_definitions(
  maire_hydro_2d_regression
  PRIVATE
  FLECSI_SP_BURTON_MESH_DIMENSION=2
  FLECSI_SP_BURTON_MESH_EXTRAS
  FLECSI_ENABLE_SPECIALIZATION_TLT_INIT
  FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
)

create_regression_test(
  NAME flecsale_sedov_2d_1procs
  COMMAND
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:maire_hydro_2d_regression> ${MPIEXEC_POSTFLAGS} -m
    ${FLECSALE_DATA_DIR}/meshes/sedov_32x32.g
    -ll:csize 9000
  COMPARE
    sedov_2d-solution_rank000000.txt
  STANDARD
    ${CMAKE_CURRENT_SOURCE_DIR}/sedov_2d-solution_rank000000.txt.1procs.std
)
