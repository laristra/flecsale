# ALE - NGC Subproject

# Summary
This is a *hopefully* parallel unstructured solver infrastructure for devlopping multi-physics simulations. The plan is to support 2D and 3D arbitrary polyhedral meshes distributed over hundreds to thousands of nodes. The code is built upon [FleCSI](https://github.com/flecsi/flecsi), an open source general purpose set of tools for execution and state control.

# Thirdparty Libraries
This project uses a number of third party libraries.  Below is a list of the required thirdparty libraries and their respective dependants. The library dependency graph is as follows:
* [Cinch](https://github.com/losalamos/cinch) - A set of utilities and configuration options to simplify [CMake](https://cmake.org/) configuration.
* [FleCSI](https://github.com/flecsi/flecsi) - For mesh/state and execution control.
  - [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) - A graph partitioner
  - [Scotch](https://www.labri.fr/perso/pelegrin/scotch/) - A graph partitioner
  - [ExodusII](https://sourceforge.net/projects/exodusii/) - A finite-element data file format
    - [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) - Provides machine-independant file formats
      - [HDF5](https://www.hdfgroup.org/HDF5/) - A data model and file format

# Installation
This project uses [Git](https://git-scm.com/) for revision control and distribution, and [CMake](https://cmake.org/) for build configuration.  Below are some general instructions for obtaining and building ALE.

First, make sure that [CMake](https://cmake.org/) and [Git](https://git-scm.com/) are available on your system.  Typing
 
    which cmake git 
   
should correctly locate the two programs and return something like

    /usr/bin/cmake
    /usr/bin/git

If it does not, and you get something like *"which: no gits in (...)"*, try loading the necessary modules.  On most machines, the modules are loaded using

    module load cmake git

To clone the repository, type

    git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/ngc/ale.git
    
**Make sure to include the *"\-\-recursive"* so that all of the submodules are cloned as well.**  Once the code is checked out, create a build directory somewhere. 

    mkdir build
    cd build
    
Building the code is a three step process:  1) obtaining the thirdparty libraries; 2) building the thirdparty libraries; and 3) compiling the final ALE project.

## 1. Obtain the Thirdparty Libraries
If you already have the necessary libraries in a folder somewhere, you can skip this step.  Otherwise, the build system can download them for you.  A good place to put them is inside the reposotory itself.  An empty folder was created inside the git repository for this purpose.  It is located in **<ALE_DIR>/thirdparty/files** where **<ALE_DIR>** is the location of the cloned git repository.

    TPL_DOWNLOAD_PATH=<ALE_DIR>/thirdparty/files <ALE_DIR>/arch/download-tpl.sh

## 2. Build the Thirdparty Libraries
To build the thirdparty libraries, you need a C and C++ compiler installed.  You might have to load the necessary modules

    module load gcc/5.3.0

sdf

    mkdir tpl
    cd tpl
    /path/to/ale/

asdfasdf
