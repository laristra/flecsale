set(SCOTCH_NAME scotch)

if (NOT ENABLE_SCOTCH)
  return()
endif()

find_program(FLEX flex)
if (NOT FLEX)
   message( FATAL_ERROR
            "'flex' lexical parser not found. Cannot build scotch." )
endif()

if (APPLE)
    set(SCOTCH_MAKE Makefile.inc.i686_mac_darwin8)
    set(SCOTCH_LDFLAGS "")
    set(SCOTCH_CFLAGS "-Drestrict=__restrict -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_TIMING_OLD -DSCOTCH_PTHREAD -DSCOTCH_RENAME -DCOMMON_PTHREAD_BARRIER")
else ()
    if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        set(SCOTCH_MAKE Makefile.inc.x86-64_pc_linux2)
        set(SCOTCH_CFLAGS "-DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict -DIDXSIZE64")
    else ()
        set(SCOTCH_MAKE Makefile.inc.i686_pc_linux2)
        set(SCOTCH_CFLAGS "-DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict")
    endif ()
    set(SCOTCH_LDFLAGS "${ZLIB_LIBRARIES} -lm -lrt -lpthread")
endif ()

tpl_dowload_path( URL ${SCOTCH_URL} ${TPL_DOWNLOAD_PATH} )

ExternalProject_Add( ${SCOTCH_NAME}
  DEPENDS ${ZLIB_NAME}

 URL ${URL}/${SCOTCH_TGZ}
 URL_MD5 ${SCOTCH_MD5}
 UPDATE_COMMAND ""

 BUILD_IN_SOURCE 1

 CONFIGURE_COMMAND rm -f src/Makefile.inc

 COMMAND ln -s Make.inc/${SCOTCH_MAKE}
               src/Makefile.inc

 BUILD_COMMAND 
   $(MAKE) -C src
     "CCS=${CMAKE_C_COMPILER}"
     "CCD=${CMAKE_C_COMPILER}"
     "CFLAGS=${CMAKE_C_FLAGS} ${SCOTCH_CFLAGS}"
     "LDFLAGS=${SCOTCH_LDFLAGS}"
     "CLIBFLAGS=-fPIC"     
     scotch
 INSTALL_COMMAND 
   $(MAKE) -C src 
     "MKDIR=mkdir -p"
     prefix=${CMAKE_INSTALL_PREFIX} 
     install
)
