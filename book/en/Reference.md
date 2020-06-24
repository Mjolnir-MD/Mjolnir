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

## Including different files

Everywhere in the input file, a value corresponding to the key, `include`, has a special meaning.

A value corresponds to `include` is a string or an array of strings, specifying filenames.
The path would be specified via `input.path` value in the `[files]` table.

If `include` is specified, the values defined in those specified toml files will be merged under the table that has `include`.

In the following case,

```toml
# main.toml
[[forcefields]]
include = [
    "bond-length.toml",
    # other forcefields ...
]
```

```toml
# bond-length.toml
[[local]]
interaction = "BondLength"
potential   = "Harmonic"
parameters  = [
    # ...
]
```

The file internally becomes to the following.

```toml
[[forcefields]]
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
parameters  = [
    # ...
]
```

Note that Mjolnir expands those include files only once.
