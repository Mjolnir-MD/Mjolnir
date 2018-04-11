Mjolnir
==========

[![Build Status](https://travis-ci.org/ToruNiina/Mjolnir.svg?branch=master)](https://travis-ci.org/ToruNiina/Mjolnir)
![spec](https://img.shields.io/badge/spec-unstable-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Molecular Dynamics Simulation Software written in c++11.

## Description

Flexible, well-organized, and modern molecular dynamics simulation code.

Mainly focused on Coarse-Grained MD simulation.

### Goals

1. Flexibility: To make it easy to implement new forcefields.
2. Reliability: To do everything that a user wants and nothing more than that.
3. Efficiency: To be fast whenever possible without breaking the above goals.

## Build

Testing code depends on the Boost.unittest framework.
If you want to run test codes, please install boost.

```sh
$ git submodule init   # (if you've just cloned this repo)
$ git submodule update # (ditto)
$ mkdir build          # (ditto)
$ cd build
$ cmake ..
$ make
$ make test # optional
```

`g++-5` and `clang++-3.7` on `Ubuntu 14.04`, `clang++` on `OS X` are tested on __Travis CI__.

## Example

__NOTE__: Specifications of this project is currently unstable.
Most of the following procedures may be changed.

Example input files are in the `input` directory.
You can run them with the command below.

```console
$ mjolnir input/sh3_AICG2+.toml
```

## Author

* Toru Niina
  * the original designer and implementer.
* Yutaka Murata
  * adds implicit membrane potential.

## Licensing terms

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2016-2018 Toru Niina

All rights reserved.
