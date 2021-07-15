+++
title = "ForceFields"
weight = 600
bookCollapseSection = true
+++

# `[[forcefields]]`

`ForceField` is a set of potential energy functions applied to the particles.

The units of the parameters are defined in [`[units]`]({{<relref "/docs/reference/units">}}).

To support simulation methods that uses several forcefields simultaneously, `[[forcefields]]` is defined as an array of tables.
But normally, we need only one forcefield per simulation.

## env

In `[[forcefields]]` tables, a special field, `env`, can be used.

By defining a variable under the `env` table, you can name the value.
Named values can be referenced via its name in the parameter section.

```toml
[[forcefields.local]]
env.pi = 3.1416 # this substitutes `3.1416` to `"pi"`
parameters = [
    {indices = [0, 1], k = 100.0, v0 = "pi"},
    {indices = [1, 2], k = 100.0, v0 = "pi"},
    {indices = [2, 3], k = 100.0, v0 = "pi"},
]
```

This feature can be used with [file inclusion feature]({{<relref "/docs/reference/_index.md#including-different-files">}}).

## [LocalForceFiled]({{<relref "local">}})

A set of interactions that is applied to a specific set of particles.
E.g. bond length, bond angle, or dihedral angle interaction.

## [GlobalForceFiled]({{<relref "global">}})

A set of interactions that is applied to all the possible pair of particles.
E.g. electrostatic, or Lennard-Jones interaction.

## [ExternalForceFiled]({{<relref "external">}})

A set of interactions that is applied to particles from external environment.
E.g. restraint on a specific position, or a box interaction.

## [ConstraintForceField]({{<relref "constraint">}})

A set of constraints. Since it affects on the topology, it is defined as a forcefield.

## [MultipleBasinForceField]({{<relref "MultipleBasinForceField.md">}})

A special meta-forcefield that allows to concatenate two or more forcefields smoothly.

## [HybridForceField]({{<relref "HybridForceField.md">}})

A meta-forcefield that is a linear combination of two different forcefields.
