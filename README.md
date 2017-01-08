Mjolnir
====

![spec](https://img.shields.io/badge/spec-unstable-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Molecular Dynamics Simulation Software written in c++11.

## build

Mjolnir depends on [TOMLParser](https://github.com/ToruNiina/TOMLParser).
Before build, please Download the library and add include path.

```sh
$ mkdir build
$ cd build
$ cmake ..
$ make
```

## testing

After build, you can run test codes. 

```sh
$ make test
```
