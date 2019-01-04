# Mjolnir

[![Build Status](https://travis-ci.org/ToruNiina/Mjolnir.svg?branch=master)](https://travis-ci.org/ToruNiina/Mjolnir)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Molecular Dynamics Simulation Software written in C++11.

## Description

Flexible, well-organized, and modern molecular dynamics simulation code.

Mainly focused on Coarse-Grained MD simulation.

### Goals

 1. Flexibility: To make it easy to implement new forcefields.
 2. Reliability: To do everything that a user wants and nothing more than that.
 3. Efficiency: To be fast whenever possible without breaking the above goals.

## Build

All the depending libraries are automatically downloaded in the CMake script.

```console
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test # optional
```

After this, you will find executable binary in `bin` directory.

The codes are tested with the following compilers on [Travis CI](https://travis-ci.org/ToruNiina/Mjolnir).
  - `g++-7`, `g++-6`, `g++-5` on `Linux`
  - `clang++-7`, `clang++-6`, `clang++-5` on `Linux`
  - `clang++` on `OS X`

## Dependency

Mjolnir depends on [toml11](https://github.com/ToruNiina/toml11).
The test codes depend on [Boost C++ Library](https://www.boost.org/).

These libraries are automatically downloaded in the CMake script.
You need nothing to install.

If you have already installed recent version of Boost Library (1.67.0 or later),
you can use it by passing `-DBOOST_ROOT=/path/to/boost` option to `cmake`.
It will save your time.

## Example

Example input files are in the `input` directory.
You can run them with the command below.

```console
$ ./bin/mjolnir input/sh3_AICG2+.toml
```

## Author

  - Toru Niina
    - the original designer and implementer.
  - Yutaka Murata
    - adds implicit membrane potential.

## Licensing terms

This product is licensed under the terms of the [MIT License](LICENSE).

  - Copyright (c) 2016-2018 Toru Niina

All rights reserved.
