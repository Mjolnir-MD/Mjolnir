Mjolnir
==========

[![Build Status](https://travis-ci.org/ToruNiina/Mjolnir.svg?branch=master)](https://travis-ci.org/ToruNiina/Mjolnir)
![spec](https://img.shields.io/badge/spec-unstable-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Molecular Dynamics Simulation Software written in c++11.

Build
----------

__Mjolnir__ depends on [TOMLParser](https://github.com/ToruNiina/TOMLParser).

```sh
$ git submodule init   # (if you've just cloned this repo)
$ git submodule update # (ditto)
$ mkdir build          # (ditto)
$ cd build
$ cmake ..
$ make
```

Testing
----------

After build, you can run test codes.

```sh
$ make test
```

## Author

* Toru Niina

## Licensing terms

This product is licensed under the terms of the [MIT License](LICENSE).

- Copyright (c) 2016- Toru Niina

All rights reserved.
