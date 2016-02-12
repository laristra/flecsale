# ALE - NGC Subproject

# Summary
This is a *hopefully* parallel unstructured solver infrastructure for
devlopping multi-physics simulations. The plan is to support 2D and 3D
arbitrary polyhedral meshes distributed over hundreds to thousands of
nodes. The code is built upon
[FleCSI](https://github.com/flecsi/flecsi), an open source general
purpose set of tools for execution and state control.


# Thirdparty Libraries
This project uses a number of third party libraries.  Below is a list
of the required thirdparty libraries and their respective
dependants. The library dependency graph is as follows:

- [Cinch](https://github.com/losalamos/cinch) - A set of utilities and
  configuration options to simplify [CMake](https://cmake.org/)
  configuration.
  
  
- [FleCSI](https://github.com/flecsi/flecsi) - For mesh/state and
  execution control. 
    - [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) -
      A graph partitioner
      
    - [Scotch](https://www.labri.fr/perso/pelegrin/scotch/) - A graph
      partitioner
      
    - [ExodusII](https://sourceforge.net/projects/exodusii/) - A
      finite-element data file format
        - [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) -
          Provides machine-independant file formats
            - [HDF5](https://www.hdfgroup.org/HDF5/) - A data model and
              file format

# Installation
This project uses [Git](https://git-scm.com/) for revision control and
distribution, and [CMake](https://cmake.org/) for build configuration.
Below are some general instructions for obtaining and building ALE.

First, make sure that [CMake](https://cmake.org/) and
[Git](https://git-scm.com/) are available on your system.  Typing
 
    which cmake git 
   
should correctly locate the two programs and return something like

    /usr/bin/cmake
    /usr/bin/git

If it does not, and you get something like *"which: no gits in
(...)"*, try loading the necessary modules.  On most machines, the
modules are loaded using

    module load cmake git

To clone the repository, type

    git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/ngc/ale.git
    
**Make sure to include the *"\-\-recursive"* so that all of the
submodules are cloned as well.**  Once the code is checked out,
create a build directory somewhere for the thirdparty libraries and
the ALE code

    mkdir -p build/tpl
    mkdir -p build/ale
    ls build
    $ ale tpl
    
Building the code is a two step process:  1) building the thirdparty
libraries and 2) compiling the final ALE project.

## 1. Build the Thirdparty Libraries
The thirdparty libraries have already been downloaded for you on most
machines.  To build the thirdparty libraries, you need a C and C++
compiler installed.  You might have to load the necessary modules

    module load gcc/5.3.0

Make a directory for the thirdparty libraries in **build** and compile
the libraries.  The following commad will build and install them in a
folder **build/tpl/install**:

    cd build/tpl
    CC=gcc CXX=g++ TPL_INSTALL_PREFIX=./install <ALE_DIR>/arch/cmake-tpl.sh
    make -j

## 2. Build the ALE Project

Now build the final ALE project.  Go into the **build/ale** directory
we created, run cmake, then make the code.  Execute the following:

    cd ../ale
    CC=gcc CXX=g++ TPL_INSTALL_PREFIX=../tpl/install <ALE_DIR>/arch/cmake-ale.sh
    make -j
    
You do not have to rebuild the thirdparty libraries if you make
changes to the ALE source code in **\<ALE_DIR\>**.  Simply recompile the
ALE source using the following:

    cd build/ale
    make -j

# Code Structure

