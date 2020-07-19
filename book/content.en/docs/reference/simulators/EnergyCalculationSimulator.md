+++
title  = "EnergyCalculation"
weight = 5000
+++

# EnergyCalculation

It does not perform any simulation.

It calculates energy of each snapshot written in a given trajectory file, taking one [`ForceField`]({{<relref "/docs/reference/forcefields">}}).

It does not use any [`Integrator`]({{<relref "/docs/reference/integrators">}}).
Also, coordinates specified in [`[[systems]]`]({{<relref "/docs/reference/system">}}) will be overwritten by the trajectory file.

Still, it requires [`[[systems]]`]({{<relref "/docs/reference/system">}}) field because group of the particle is given in `[[systems]]`.

## Example

```toml
[simulator]
type          = "EnergyCalculation"
boundary_type = "PeriodicCuboid"
precision     = "double"
parallelism   = "OpenMP" # optional
file          = "example_position.dcd"
```

## Input Reference

- `type`: String
  - Name of the simulator. Here, it is `"EnergyCalculation"`.
- `boundary_type`: String
  - Type of the boundary condition. The size will be specified in [`[[systems]]`]({{<relref "/docs/reference/system">}}).
  - `"Unlimited"`: No boundary condition will applied.
  - `"Periodic"`: Periodic boundary condition will be applied. The shape is recutangular box.
- `precision`: String
  - Precision of floating point number used in the simulation.
  - `"float"`: 32bit floating point number.
  - `"double"`: 64bit floating point number.
- `parallelism`: String (Optional. By default, `"sequencial"`.)
  - `"OpenMP"`: OpenMP implementation will be used.
  - `"sequencial"`: It runs on single core.
- `file`: String
  - Trajectory file. It considers the input path specified in [`[files]`]({{<relref "/docs/reference/files">}}).
