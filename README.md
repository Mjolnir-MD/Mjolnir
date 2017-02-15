Mjolnir
==========

[![Build Status](https://travis-ci.org/ToruNiina/Mjolnir.svg?branch=master)](https://travis-ci.org/ToruNiina/Mjolnir)
![spec](https://img.shields.io/badge/spec-unstable-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Molecular Dynamics Simulation Software written in c++11.

## Purpose

Flexible, well-organized, and modern molecular dynamics simulation code.

### Goals

1. __Mjolnir__ should be able to execute arbitrary molecular model and arbitrary force field.
  - in other words, it should be easy to add a novel molecular model or force field to __Mjolnir__.
2. __Mjolnir__ should do everything that user requires and should not do anything more than that.
3. Whenever possible without breaking the above goals, __Mjolnir__ should be fast.

## Build

Testing code depends on Boost.unittest framework. please install boost.unittest.

`g++-5` and `clang++-3.7` on `Ubuntu 14.04 Trusty Tahr` are tested on __Travis CI__.

```sh
$ git submodule init   # (if you've just cloned this repo)
$ git submodule update # (ditto)
$ mkdir build          # (ditto)
$ cd build
$ cmake ..
$ make
```

## Example

__NOTE__:
Specification of this project is currently unstable.
Most of the following procedures (mainly unhandy or confusing steps) may be changed.

### SRC tyrosine kinase SH3 domain

download `1SRL.pdb` from Protein Data Bank.

First, decide some parameters and write them to `toml` file.

```toml
# 1SRL.toml file
[simulator]
file_name        = "1SRL"                 # name of output files
time_integration = "Underdamped Langevin" # time integration method
boundary         = "Unlimited"            # boundary condition
delta_t          = 0.4                    # delta t
total_step       = 50_000                 # number of total steps
save_step        = 100                    # number of time interval to output
temperature      = 300.0                  # it is just a temperature.
kB               = 1.986231313e-3         # value of boltzmann constant
seed             = 2374                   # seed of random variable.
```

__Jarngreipr__ on current version generates Clementi-Go forcefield in default.

```sh
$ jarngreipr 1SRL.pdb >> 1SRL.toml
$ mjolnir 1SRL.toml
```

simulation takes about 2~3 seconds. after that, you can see the trajectory
using molecular visualization software like `vmd`, `pymol` or others.

```sh
$ [visualizer] 1SRL.xyz
```

### Beta-Galactosidase in complex with allolactose

download `1JZ8.pdb` from Protein Data Bank.

Then extract `ATOM` line in pdb file.

```sh
$ grep "^ATOM" 1JZ8.pdb > 1JZ8_ATOM.pdb
```

After that, the procedure is same as the above one.
You can run simulation of protein complex with ease.

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
