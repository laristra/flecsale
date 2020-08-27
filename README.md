[![Build Status](https://travis-ci.com/laristra/flecsale.svg?branch=master)](https://travis-ci.com/laristra/flecsale)
[![pipeline status](https://gitlab.lanl.gov/laristra/flecsale/badges/master/pipeline.svg)](https://gitlab.lanl.gov/laristra/flecsale/pipelines)

# FleCSALE

FleCSALE is a computer software package developed for studying problems that
can be characterized using continuum dynamics, such as fluid flow. It is
specifically developed for existing and emerging large distributed memory
system architectures. FleCSALE uses the Flexible Computational Science
Infrastructure ([FleCSI](https://github.com/laristra/flecsi)) project for mesh
and data structure support. FleCSALE has support for multi-phase fluids and
tabular equation of state (EOS).

# Requirements

## Minimal

- [Boost C++ Libraries](boost.org) >= 1.56
- C++17 compliant compiler  (gcc >= 7.0, clang>=3.9)
- [CMake](http://www.cmake.org/) >= 3.0
- [Exodus](https://github.com/gsjaardema/seacas) to read/write ExodusII formatted files
- [FleCSI](https://github.com/laristra/flecsi)
- MPI implementation
- [ParMETIS](https://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) for mesh partitioning

## Optional

- [Doxygen](http://doxygen.org) to generate documentation
- [Lua](http://lua.org) >= 5.2 for application input
- [OpenSSL](http://openssl.org) for field data checksums
- [Python](http://python.org) >= 2.7 for regression tests
- [VTK](http://vtk.org) to read/write VTK formatted files


# Getting the code

This project uses [Git](https://git-scm.com/) for revision control and
distribution, and [CMake](https://cmake.org/) for build configuration.
Below are some general instructions for obtaining and building FleCSALE.

FleCSALE uses git submodules, so it must be checked out recursively.  Type

    $ git clone --recursive git@github.com:laristra/flecsale.git
    
to clone the repository using ssh, or 

    $ git clone --recursive https://github.com/laristra/flecsale.git
    
to clone using https.
    
**Make sure to include the `--recursive` so that all of the
submodules are cloned as well.** 

# Getting the dependencies

If you are LANL internal users on LANL machines, please refer to the README on mm.

Make sure you have at least the minimum compiler version in your environment.

Now, you need to download Spack if you don't already have one.
```
$ cd ${PROJECT_DIR}
$ git clone https://github.com/spack/spack.git
```
Then do
```
$ source ${PROJECT_DIR}/spack/share/spack/setup-env.sh
$ spack compiler find --scope site
==> Added 1 new compilers to ${PROJECT_DIR}/spack/etc/spack/compilers.yaml
    gcc@<compiler-version>
==> Compilers are defined in the following files:
    ${PROJECT_DIR}/spack/etc/spack/compilers.yaml

$ spack compiler list
==> Available compilers
-- gcc <os>-<target> -------------------------------------------
gcc@<compiler-version>
```
to get Spack into your environment and see what compilers you have that
Spack can find automatically. 

Next, we need to get `ristra_spackages` and add the folder that contains custom ristra spackages to Spack
```
$ git clone https://github.com/laristra/ristra_spackages.git
$ spack repo add --scope site ristra_spackages/spack-repo
==> Added repo with namespace 'lanl_ristra'.

$ spack repo list
==> 2 package repositories.
lanl_ristra        ${PROJECT_DIR}/ristra_spackages/spack-repo
builtin            ${PROJECT_DIR}/spack/var/spack/repos/builtin
```
Now, assuming you have the compiler you want recognized by Spack
and added the custom spack-repo, you could just do the install for a mpi backend
using mpich like this
```
$ mkdir spack_env
$ spack env create --without-view -d spack_env
$ spack -e spack_env install --show-log-on-error flecsalemm-deps+hdf5%gcc@<compiler-version> backend=mpi ^mpich@3.2.1+slurm
```
to get all the dependencies and all their dependencies installed from
scratch.

But if you want to save time or there is some packages that spack
has trouble installing, you could let Spack know what
packages or modules you already have on the system by adding
`packages.yaml` to your `${PROJECT_DIR}/spack/etc/spack/linux`, which could look something like this:
```
packages:
  perl:
    paths:
      perl@5.16.3: /usr
  numactl:
    paths:
      numactl@system: /usr
  python:
    modules:
      python@2.7.3: python/2.7.3
      python@3.5.1: python/3.5.1
  mpich:
    modules:
      mpich@3.2.1%gcc@7.3.0: mpich/3.2.1-gcc_7.3.0
  openmpi:
    modules:
      openmpi@3.1.4%gcc@7.3.0: openmpi/3.1.4-gcc_7.3.0
  cmake:
    modules:
      cmake@3.12.4: cmake/3.12.4
```
Then the installation from Spack will take less time.

After Spack finishes the installation, you can load them into
your environment by doing
```
$ spack -e spack_env module tcl refresh -y
$ spack -e spack_env env loads -r
$ source spack_env/loads
```

# Building

Building the code is simply performed through the following steps
using [CMake](https://cmake.org/):

    $ mkdir build
    $ cd build
    $ CC=gcc CXX=g++ cmake /path/to/source/directory -DFLECSI_RUNTIME_MODEL="mpi" -DENABLE_COLORING=on [options]
    $ make -j

The environment variables `CC` and `CXX` are only necessary to select
specific compilers and may be omitted.  The first parameter provided
to [CMake](https://cmake.org/) must be the root of the source
directory created by cloning the FleCSALE repository.

Options provded to the [CMake](https://cmake.org/) command line can be
any CMake build options listed below (use `-Doption_name=value` to
specify an option.  For example, to build the unit tests, specify
`-DENABLE_UNIT_TESTS=ON`.

# CMake installation options


 - `CMAKE_BUILD_TYPE`:  Type of build: `Release` (for users) or `Debug` (for developers)
 - `ENABLE_DOXYGEN`:  Generate HTML API documentation with Doxygen - Default is `OFF`
 - `ENABLE_LUA`: Enable application input with Lua - Defaults to `ON` if Lua was found
 - `ENABLE_REGRESSION_TESTS`: Build the regression tests - Defaults to `ENABLE_UNIT_TESTS`
 - `ENABLE_UNIT_TESTS`:  Build the unit tests - Default is `OFF`
 - `FLECSI_RUNTIME_MODEL`: Parallel backend: `mpi` (for most users), `legion`, or `hpx`

# Release

This software has been approved for open source release and has
been assigned **LA-CC-16-076**.

# Copyright

Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:  

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
