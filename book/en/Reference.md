# Reference

Mjolnir does not have any command-line option. It uses a toml file as the input.
So a command to execute simulation looks like the following.


```sh
$ ./bin/mjolnir sample.toml
```

This document introduces how to write an input file. Here, a brief introduction
of major components in an input file are provided.

For detailed information about TOML v0.5.0, the file format used in Mjolnir,
can be found in the [toml-lang official repository](https://github.com/toml-lang/toml).

An input file for Mjolnir has of 5 main components.
Each of them is a table or an array of tables.

## [`[files]`](refs/Files.md)

It specifies the output path, filename, and format.
It also specify input path if required.

## [`[units]`](refs/Units.md)

It specifies a set of units to be used, e.g. `angstrom`, `kJ/mol`, etc.

The values of energies, positions, distances that would be provided later,
are considered to be expressed in this units.

## [`[simulator]`](refs/Simulator.md)

It specifies a protocol of a simulation.

E.g., which floating-point would be used, what boundary condition would be used,
how many timesteps does it simulate, etc.

You can provide this table as another toml file by specifying only a filename.

## [`[[systems]]`](refs/System.md)

It specifies a set of particles, boundary shape, system parameters such as
temperature, ionic strength, etc.

It is defined as an array of table for some kinds of simulations that uses
several replicas of simulation box. But in the normal case, only one
`[[systems]]` are required.

You can provide this table as another toml file by specifying only a filename.

## [`[[forcefields]]`](refs/ForceField.md)

It specifies forcefield parameters.

It is defined as an array of table for some kinds of simulations that uses
several different forcefields. But in the normal case, only one `[[forcefields]]`
are required.

You can provide this table as another toml file by specifying only a filename.
