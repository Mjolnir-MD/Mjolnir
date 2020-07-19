+++
title = "Constraint"
weight = 4000
+++

# ConstraintForceField

`ConstraintForceField` is a set of consraints on a distance between two particles.

Currently, only [`g-BAOABLangevin`]({{<relref "/docs/reference/integrators/gBAOABLangevinIntegrator.md">}}) integrator can handle this.

## Example

```toml
[[forcefields.constraint]]
topology      = "bond"
max_iteration = 500
tolerance     = 1e-6
parameters  = [
    {indices = [0, 1], offset = 100, v0 = 3.8},
    # ...
]
```

## Input Reference

- `topology`: String
  - Name of this constraint on [`"Topology"`]({{<relref "/docs/reference/forcefields/Topology.md">}})
- `max_iteration`: Integer
  - It repeats the correction up to this number of times until convergence.
  - If it does not converge, a warning will be shown. In most cases, it will diverge.
- `tolerance`: Floating
  - Tolerance to be used in convergence check.
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 2)
    - The index of particles to be constrained. The index is 0-based.
  - `offset`: Integer(Optional. By default, 0.)
    - Offsets of index.
  - `v0`: Floating
    - The reference distance between particles.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
