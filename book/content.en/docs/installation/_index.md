+++
title  = "Installation"
weight = 10
+++

# Installation

This page describes how to build Mjolnir.

## Prerequisities

- linux or Unix (e.g. OS X)
- C++11 compatible compiler
- git
- Make
- CMake
- Boost C++ Library

### Operating System

It assumes that the file pathes conform to the posix standard.

Also, automatic tests only run on Linux and OS X. Thus it is not guaranteed
to work on Windows properly.

### C++ compiler

It requires C++11 compatible compiler. Since now it is 2021, default compilers
in most of the systems are compatible to C++11.

But I recommend to use the later versions of compilers. Compilers are also
in development, so the later versions generally have less bugs and generate
faster executables.

### Git

Mjolnir contains `toml11` library as a `git-submodule`.
 Git is required to download it.

### CMake

Mjolnir uses CMake as a build system.

The `CMakeLists.txt`s expects relatively later versions of CMake.
So it is recommended to use the later versions of CMake.

### Boost

Test codes of Mjolnir depends on the Boost library.

Although `CMakeLists.txt` automatically downloads Boost, it makes build process
faster to use a library that is already built.

## Building

To build Mjolnir, run the following commands.

```console
$ git clone https://github.com/Mjolnir-MD/Mjolnir.git
$ cd Mjolnir
$ mkdir build
$ cd build
$ cmake .. # commonly-used options are listed below.
$ make
$ make test
```

### Options for CMake

All those are optional variables.

- `-DCMAKE_CXX_COMPILER=/path/to/compiler`
  - It is an option for CMake. You can choose what compiler would be used.
- `-DUSE_OPENMP=(ON|OFF)`
  - `ON` by default.
  - If `ON`, compile with OpenMP (if it is available).
    Then you can parallelize your simulation with OpenMP.
- `-DFIND_BOOST=(ON|OFF)`
  - `OFF` by default.
  - If `ON`, CMake looks boost library that is already installed.
    If it does not exist, the build fails.
- `-DBOOST_ROOT=/path/to/boost`
  - It is an option for `FindBoost` package in CMake.
    You can specify the path to boost you want to use.
- `-DSEPARATE_BUILD=(ON|OFF)`
  - `OFF` by default.
  - If `ON`, compile codes separately and link them after compilation.
    This feature is for developers who compiles it a lot of times.
- `-DSINGLE_PRECISION_SUPPORT=(ON|OFF)`
  - `ON` by default.
  - If `OFF`, `float` version of the code will be ignored and the compilation finishes quickly.
- `-DDOUBLE_PRECISION_SUPPORT=(ON|OFF)`
  - `ON` by default.
  - If `OFF`, `double` version of the code will be ignored and the compilation finishes quickly.
- `-DUNLIMITED_BOUNDARY_SUPPORT=(ON|OFF)`
  - `ON` by default.
  - If `OFF`, `UnlimitedBoundary` would not be able to be used, but the compilation finishes quickly.
- `-DPERIODIC_BOUNDARY_SUPPORT=(ON|OFF)`
  - `ON` by default.
  - If `OFF`, `CuboidalPeriodicBoundary` would not be able to be used, but the compilation finishes quickly.

