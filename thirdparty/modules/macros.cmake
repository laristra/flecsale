################################################################################
# select to download or use existing
################################################################################
function(tpl_dowload_path  path package_url )
  if (ARGV2) 
    set( ${path} ${ARGV2} PARENT_SCOPE )
  else()
    set( ${path} ${package_url} PARENT_SCOPE )
  endif()
endfunction()




################################################################################
# Download libraries
################################################################################
function(download_libraries  path )

  include(versions)

  set( packages 
    szip 
    zlib 
    hdf5 
    netcdf 
    exodus 
    metis
    scotch )

  foreach ( pkg ${packages} )
    message( STATUS "Downloading ${pkg}" )

    # get the file names
    string(TOUPPER ${pkg} PKG)
    set( filename ${${PKG}_TGZ} )
    set( url "${${PKG}_URL}/${filename}" )
    set( tgt_file ${path}/${filename} )
    set( md5 ${${PKG}_MD5} )

    set( download TRUE )
    
    # check if it exists and the checksum is correct
    if( EXISTS "${tgt_file}" ) 
      message( STATUS "${pkg} exists, verify checsum." )
      # generate new checksum
      execute_process( 
        COMMAND ${CMAKE_COMMAND} -E md5sum ${tgt_file}
        OUTPUT_VARIABLE CalculatedCheckSum 
        OUTPUT_STRIP_TRAILING_WHITESPACE )
      separate_arguments( CalculatedCheckSum )      
      list(GET CalculatedCheckSum 0 CalculatedMD5)
      # check the checksum
      if ("${md5}" STREQUAL "${CalculatedMD5}")
        message( STATUS "checsum matches, nothing to do." )
        set( download FALSE )
      else()
        message( STATUS "Checsum mismatch; expected = ${md5}, calculated = ${CalculatedMD5}" )
      endif()
    endif()
    
    # download if needed
    if ( download ) 
      execute_process( COMMAND wget ${url} -O ${tgt_file} )
    endif()

  endforeach()


endfunction()

