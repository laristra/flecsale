#!/bin/bash
#~----------------------------------------------------------------------------~#
# placeholder
#~----------------------------------------------------------------------------~#

compiler=$1

TPL_DEFAULT_PATH=/usr/projects/ngc/public/ale-thirdparty
TPL_DEFAULT_DOWNLOAD_PATH=${TPL_DEFAULT_PATH}/src
TPL_DEFAULT_INSTALL_PREFIX=${TPL_DEFAULT_PATH}/builds/gcc/install

if [ "$compiler" == "gcc" ]; then
  echo "Using GNU environment"
  export CC=gcc
  export CXX=g++  
elif [ "$compiler" == "clang" ]; then
  echo "Using clang environment"  
  export CC=clang
  export CXX=clang++ 
elif [ "$compiler" == "own" ]; then
  echo "Using user environment"
else
  echo "no build environment specified"
  echo "Usage: $(basename $0) env_type"
  echo "env_type : gcc   - use GNU compilers"
  echo "           clang - use clang compilers"
  echo "           own   - let user set their own environment"
  exit 1
fi


#~---------------------------------------------------------------------------~-#
# placeholder
#~---------------------------------------------------------------------------~-#
