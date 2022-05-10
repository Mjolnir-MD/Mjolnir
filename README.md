# Mjolnir

[![build](https://github.com/Mjolnir-MD/Mjolnir/workflows/build/badge.svg)](https://github.com/Mjolnir-MD/Mjolnir/actions)
[![document](https://github.com/Mjolnir-MD/Mjolnir/workflows/document/badge.svg)](https://github.com/Mjolnir-MD/Mjolnir/actions)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b55282103ca74dd5b9b0022a3af99f3b)](https://www.codacy.com/app/ToruNiina/Mjolnir?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Mjolnir-MD/Mjolnir&amp;utm_campaign=Badge_Grade)
[![Latest Version](https://img.shields.io/github/release/Mjolnir-MD/Mjolnir.svg)](https://github.com/Mjolnir-MD/Mjolnir/releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

## Description

Flexible, well-organized, and modern molecular dynamics simulation package, mainly for coarse-grained models.

The detailed information can be found at https://mjolnir-md.github.io/Mjolnir/ .

### Goals

1. Flexibility: To make it easy to implement new forcefields.
2. Transparency: To do everything that a user wants and nothing more than that.
3. Efficiency: To be fast whenever possible without breaking the above goals.

### Capability

- [x] Modern NVT/NPT Langevin integrators, such as BAOAB and G-JF.
- [x] Langevin integrators with bond-length constraints based on g-BAOAB.
- [x] Several simulation algorithms, such as simulated annealing, switching forcefield, etc.
- [x] Support for all the well-known interactions (pair, bond, angle, dihedral) with many potential functions. (L-J, Debye-Hückel, harmonic, gaussian, worm-like chain, and more!)
- [x] Support for well-established Protein (AICG2+) and DNA (3SPN2/C) coarse-grained models.
- [x] Generalized multiple basin forcefield that enables conformational change in a coarse-grained model.
- [x] Parallel execution based on OpenMP.

... [and many more!](https://mjolnir-md.github.io/Mjolnir/docs/reference/)

## Build

Since Mjolnir manages depending library via git submodule, clone this repo using git and **do not** download zip or release-tarball. It will cause compilation error.

All the depending libraries are automatically downloaded in the CMake script.

```console
$ git clone https://github.com/Mjolnir-MD/Mjolnir.git
$ cd Mjolnir
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test # optional
```

After this, you will find an executable binary in `bin` directory.

The codes are automatically tested on [GitHub Actions](https://github.com/Mjolnir-MD/Mjolnir/actions).

## Dependency

Mjolnir depends on [toml11](https://github.com/ToruNiina/toml11).

The test codes depend on [Boost C++ Library](https://www.boost.org/).

These libraries are automatically downloaded in the CMake script.
You don't need to install anything manually.

If you have already installed recent version of Boost Library (1.67.0 or later), you can use it by passing `-DBOOST_ROOT=/path/to/boost` option to `cmake`.
It will save your time and system memory.

## Example

Example input files are in the `input` directory.
You can run them with the command below.

```console
$ ./bin/mjolnir input/sh3_AICG2+.toml
```

The detailed information can be found at https://mjolnir-md.github.io/Mjolnir/.

## Author

- Toru Niina
  - the original designer and implementer.
- Yutaka Murata
  - contributor implementing many features (see [CHANGELOG](CHANGELOG.md)).

## Licensing terms

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2016-2021 Toru Niina

All rights reserved.
