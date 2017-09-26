Mjolnir
==========

[![Build Status](https://travis-ci.org/ToruNiina/Mjolnir.svg?branch=master)](https://travis-ci.org/ToruNiina/Mjolnir)
![spec](https://img.shields.io/badge/spec-unstable-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Molecular Dynamics Simulation Software written in c++11.

## Purpose

Flexible, well-organized, and modern molecular dynamics simulation code.

Mainly focused on Coarse-Grained MD simulation.

### Goals

1. __Mjolnir__ should be able to execute arbitrary molecular model and arbitrary force field.
  - in other words, it should be easy to add a novel molecular model or force field to __Mjolnir__.
2. __Mjolnir__ should do everything that user requires and should not do anything more than that.
3. Whenever possible without breaking the above goals, __Mjolnir__ should be fast.

## Build

__Mjolnir__ depends on library to parse `toml` format file, [toml11](https://github.com/ToruNiina/toml11).
It is managed as a `git submodule`. Please initialize and update it after `git clone`.

Testing code depends on Boost.unittest framework. please install boost.unittest.

```sh
$ git submodule init   # (if you've just cloned this repo)
$ git submodule update # (ditto)
$ mkdir build          # (ditto)
$ cd build
$ cmake ..
$ make
```

`g++-5` and `clang++-3.7` on `Ubuntu 14.04`, `clang++` on `OS X` are tested on __Travis CI__.

## Example

__NOTE__: Specification of this project is currently unstable.
Most of the following procedures may be changed.

### SRC tyrosine kinase SH3 domain

Example input files are in `input` directory. You can run it with the command below.

```console
$ mjolnir input/sh3_AICG2+.toml
```

the trajectory will be written in `sh3_domain_unlim.xyz`.
You can see the trajectory with some protein viewer like `vmd`.

### Input File

Input file is composed of mainly 5 parts.

* general
    * define general parameter like precition, filename prefix, etc.
      for some reason, boundary condition is defined here.
* parameter
    * set physical parameters.
* system
    * define initial condition containing coordinates and system parameters
      like temperature, ionic strength, boundary, etc.
* simulator
    * setup simulator. now, only "Molecular Dynamics" is available.
* forcefield
    * define forcefields composed of pairs of interaction and potentials.

## Testing

After build, you can run test codes.

```sh
$ make test
```

Author
----------

* Toru Niina

Licensing terms
----------

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2016-2017 Toru Niina

All rights reserved.
