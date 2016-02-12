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
  module purge
  module load cmake gcc/5.3.0
  module list
  export CC=gcc
  export CXX=g++  
  TPL_DEFAULT_INSTALL_PREFIX=${TPL_DEFAULT_PATH}/builds/$compiler/install
elif [ "$compiler" == "clang" ]; then
  echo "Using clang environment"  
  module purge
  module load cmake clang/3.7.0
  module list
  export CC=clang
  export CXX=clang++ 
  export CFLAGS=--gcc-toolchain=/projects/opt/centos7/gcc/5.3.0
  export CXXFLAGS=--gcc-toolchain=/projects/opt/centos7/gcc/5.3.0
  TPL_DEFAULT_INSTALL_PREFIX=${TPL_DEFAULT_PATH}/builds/$compiler/install
elif [ "$compiler" == "own" ]; then
  echo "Using user environment"
  module list
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
