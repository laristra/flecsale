#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

return()

#------------------------------------------------------------------------------#
# Add a rule to build the executable
#------------------------------------------------------------------------------#

add_executable(simple simple.cc)

#------------------------------------------------------------------------------#
# Add link dependencies
#------------------------------------------------------------------------------#

target_link_libraries(simple ${CINCH_LINK_LIBRARIES})


#------------------------------------------------------------------------------#
# Create drivers
#------------------------------------------------------------------------------#

file(GLOB dirs *_driver)

foreach(dir ${dirs})

  if (IS_DIRECTORY ${dir})

     get_filename_component(base_name ${dir} NAME)

     file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/${base_name})
     file(COPY ${base_name} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     string(REPLACE "_" ";" temp_list ${base_name})
     list(REMOVE_AT temp_list 0)
     string(REPLACE ";" "_" target_name "${temp_list}")

     message( STATUS "Adding driver ${target_name}" )
     add_custom_target(
        ${target_name} ALL
        ${CMAKE_BINARY_DIR}/bin/flecsi ${target_name}.h
        DEPENDS ale flecsi ${CMAKE_BINARY_DIR}/bin/flecsi
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${base_name}
     )

  endif()

endforeach()

#target_link_libraries(burton ${CINCH_LINK_LIBRARIES} ${EXTRA_LIBS})


#~---------------------------------------------------------------------------~-#
# Formatting options for emacs and vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#~---------------------------------------------------------------------------~-#
