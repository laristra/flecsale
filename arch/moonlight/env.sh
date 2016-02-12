#!/bin/bash
#~----------------------------------------------------------------------------~#
# placeholder
#~----------------------------------------------------------------------------~#

# a function to print the usage
usage(){
  echo "no build environment specified"
  echo "Usage: $(basename $0) env_type"
  echo "env_type : gcc   - use GNU compilers"
  echo "           clang - use clang compilers"
  echo "           own   - let user set their own environment"
  exit 1
}


# set some variables
compiler=$1

TPL_DEFAULT_PATH=/usr/projects/ngc/public/ale-thirdparty
TPL_DEFAULT_DOWNLOAD_PATH=${TPL_DEFAULT_PATH}/src


# make sure at least one argument
if [ "$#" -eq 0 ]; then
  usage
fi

# use your own build environment
if [ "$compiler" == "own" ]; then

  echo "Using user environment"
  module list

# use a preset one
else

  echo "Using ${compiler^^} environment"
  module purge
  module load git/2.3.3 cmake/3.0.0 user_contrib boost/1.59.0
  export BOOST_INCLUDEDIR=$BOOST_INC_DIR

  if [ "$compiler" == "gcc" ]; then
    module load gcc/5.3.0
    export CC=gcc
    export CXX=g++  
    TPL_DEFAULT_INSTALL_PREFIX=${TPL_DEFAULT_PATH}/builds/$compiler/install
  elif [ "$compiler" == "clang" ]; then
    module load clang/3.7.0
    export CC=clang
    export CXX=clang++ 
    export CFLAGS=--gcc-toolchain=/usr/projects/hpcsoft/toss2/common/gcc/5.3.0
    export CXXFLAGS=--gcc-toolchain=/usr/projects/hpcsoft/toss2/common/gcc/5.3.0
    TPL_DEFAULT_INSTALL_PREFIX=${TPL_DEFAULT_PATH}/builds/$compiler/install
  else
    echo "Unknown environment selected"
    usage
  fi

  module list

fi

#~---------------------------------------------------------------------------~-#
# placeholder
#~---------------------------------------------------------------------------~-#
