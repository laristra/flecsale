#!/bin/bash

#~----------------------------------------------------------------------------~#
# placeholder
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Get the path to the project from which this script was called
#------------------------------------------------------------------------------#


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TPL_DIR=${SCRIPT_DIR}/../../thirdparty


#------------------------------------------------------------------------------#
# load environment
#------------------------------------------------------------------------------#

source ${SCRIPT_DIR}/env.sh

#------------------------------------------------------------------------------#
# Check required environment variables
#------------------------------------------------------------------------------#

if [ -z "${TPL_INSTALL_PREFIX}" ] ; then
  echo "You must set TPL_INSTALL_PREFIX in your environment"
  exit 1 
fi

if [ -z "${TPL_DOWNLOAD_PATH}" ] ; then
  #TPL_DOWNLOAD_PATH=${TPL_DIR}/files
  TPL_DOWNLOAD_PATH=${TPL_DEFAULT_DOWNLOAD_PATH}
  echo "TPL_DOWNLOAD_PATH not set, using ${TPL_DOWNLOAD_PATH}"
  echo "To change, set TPL_DOWNLOAD_PATH in your environment"
fi

#------------------------------------------------------------------------------#
# Call CMake command
#------------------------------------------------------------------------------#


cmake \
    -DCMAKE_INSTALL_PREFIX=${TPL_INSTALL_PREFIX} \
    -DTPL_DOWNLOAD_PATH=${TPL_DOWNLOAD_PATH} \
    -DTPL_DOWNLOAD_ONLY=OFF \
    -DENABLE_EXODUS=ON \
    -DENABLE_METIS=ON \
    -DENABLE_SCOTCH=ON \
    -DENABLE_TECIO=ON \
    ${TPL_DIR}

#------------------------------------------------------------------------------#
# vim: syntax=sh
#------------------------------------------------------------------------------#

#~---------------------------------------------------------------------------~-#
# placeholder
#~---------------------------------------------------------------------------~-#
