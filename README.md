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

- [FleCSI](https://github.com/laristra/flecsi)
- [cereal](https://github.com/USCiLab/cereal)
- [CMake](http://www.cmake.org/) >= 2.8
- C++14 compliant compiler  (gcc >= 5.3.0, clang>=3.7.0)

## Optional

- [Doxygen](http://doxygen.org) to generate documentation
- [Exodus](https://github.com/gsjaardema/seacas) to read/write ExodusII formatted files
- [Lua](http://lua.org) >= 5.2 for application input
- [OpenSSL](http://openssl.org) for field data checksums
- [Python](http://python.org) >= 2.7 for regression tests
- [VTK](http://vtk.org) to read/write VTK formatted files


# Getting the code

This project uses [Git](https://git-scm.com/) for revision control and
distribution, and [CMake](https://cmake.org/) for build configuration.
Below are some general instructions for obtaining and building FleCSALE.

FleCSALE uses git submodules, so it mush be checked out recursively.  Type

    $ git clone --recursive git@github.com:laristra/flecsale.git
    
to clone the repository using ssh, or 

    $ git clone --recursive https://github.com/laristra/flecsale.git
    
to clone using https.
    
**Make sure to include the `--recursive` so that all of the
submodules are cloned as well.** 

# Getting Cereal

[FleCSI](https://github.com/laristra/flecsi) uses the 
[cereal](https://github.com/USCiLab/cereal) library for serializing 
data.  To download it using git, type

    $ git clone git@github.com:USCiLab/cereal.git

to clone the repository using ssh, or 

    $ git clone https://github.com/USCiLab/cereal.git

to clone using https.

# Installation

Building the code is simply performed through the following steps
using [CMake](https://cmake.org/):

    $ mkdir build
    $ cd build
    $ CC=gcc CXX=g++ cmake /path/to/source/directory \
      -DCereal_INCLUDE_DIR=/path/to/cereal/include/ [options]
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
 - `ENABLE_OPENSSL`: Enable checksum reporting - Defaults to `ON` if OpenSSL was found
 - `ENABLE_REGRESSION_TESTS`: Build the regression tests - Defaults to `ENABLE_UNIT_TESTS`
 - `ENABLE_UNIT_TESTS`:  Build the unit tests - Default is `OFF`

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
