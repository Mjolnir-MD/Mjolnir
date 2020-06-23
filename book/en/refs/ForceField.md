# ForceField

`ForceField` is a set of potential energy functions applied to the particles.

The units of the parameters are defined in `[units]`.

To support simulation methods that uses several forcefields simultaneously,
`[[forcefields]]` is defined as an array of tables.
But normally, we need only one forcefield per simulation.

## [LocalForceFiled](LocalForceField.md)

A set of interactions that is applied to a specific set of particles.
E.g. bond length, bond angle, or dihedral angle interaction.

## [GlobalForceFiled](GlobalForceField.md)

A set of interactions that is applied to all the possible pair of particles.
E.g. electrostatic, or Lennard-Jones interaction.

## [ExternalForceFiled](ExternalForceFiled.md)

A set of interactions that is applied to particles from external environment.
E.g. restraint on a specific position, or a box interaction.

## [ConstraintForceField](ConstraintForceField.md)

A set of constraints.

## Defining `[[forcefields]]` in another file

[`include`](../Reference.md) works with `[[forcefields]]`, but there is another way to split input file, because of a historical reason.

`[[forcefields]]` tables can be defined in a different file from the main input file.

```toml
# main.toml
[files]
input.path = "./input/"

[[forcefields]]
file_name = "forcefield.toml"
```

```toml
# ./input/forcefield.toml
[[forcefields]]
# ...
```

The file path can be specified via `files.input.path`.
