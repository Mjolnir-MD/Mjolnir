# Installation

This page introduces what variables would be specified in an input file and
how it behaves.

## Prerequisities

Mjolnir requires the following stuff to build.

- linux or Unix (e.g. OS X)
- C++11 compatible compiler
- git
- Make
- CMake
- Boost C++ Library

### Operating System

Mjolnir does not use OS-specific behaviors except file paths, but it is tested
only on linux and OS X. Thus, basically, it is not guaranteed to run on any
other OS, such as Windows.

### C++ compiler

Mjolnir requires C++11 compatible compiler. Today, most of the environment
provides sufficiently new compiler by default.

If possible, it is recommended to use the latest compiler. Generally, a newer
version of compiler is smarter than the previous ones and generate the more
efficient executable binary.

### Git

Mjolnir depends on an third-party library "toml11", which is managed as a git
submodule. To download it, Git is required.

### CMake

Mjolnir uses CMake to generate a build script.

Some of the CMake modules are updated frequently. Mjolnir intends to support
newer versions of CMake rather than the older ones. So it is recommended to
update CMake to the latest version.

### Boost

Test codes of Mjolnir are written using Boost.Test library.

Mjolnir automatically downloads Boost library if it is not installed.
But pre-installed Boost library makes the building process faster.

Also, to download Boost, `wget`, `tar` and `shasum` are required. Most of the
environments have it by default.

## Building

To build Mjolnir, run the following commands.

```sh
$ git clone https://github.com/Mjolnir-MD/Mjolnir.git
$ cd Mjolnir
$ mkdir build
$ cd build
$ cmake .. # options are explained later.
$ make
$ make test
```

### Options for CMake

Useful options are listed here.

- `-DFIND_BOOST=[ON/OFF]`
  - It tries to find boost library. If it would not be found, it fails to build.
- `-DBOOST_ROOT=/path/to/boost`
  - It is an option for CMake FindBoost module. You can pass the directory where
    boost library is installed.
- `-DCMAKE_CXX_COMPILER=/path/to/compiler`
  - It is an option for CMake. You can set which compiler to use.
