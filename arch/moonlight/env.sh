#!/bin/bash
#~----------------------------------------------------------------------------~#
# placeholder
#~----------------------------------------------------------------------------~#

TPL_DEFAULT_PATH=/usr/projects/ngc/public/ale-thirdparty
TPL_DEFAULT_DOWNLOAD_PATH=${TPL_DEFAULT_PATH}/src

if [ "$1" == "gcc" ]; then
  echo "Using GNU environment"
  module purge
  module load cmake gcc/5.3.0
  module list
  CC=gcc
  CXX=g++  
elif [ "$1" == "own" ]; then
  echo "Using user environment"
  module list
else
  echo "no build environment specified"
  echo "Usage: $(basename $0) env_type"
  echo "env_type : gcc   - use GNU compilers"
  echo "           own   - let user set their own environment"
  exit 1
fi


#~---------------------------------------------------------------------------~-#
# placeholder
#~---------------------------------------------------------------------------~-#
