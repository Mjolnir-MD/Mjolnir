+++
title  = "SteepestDescent"
weight = 3000
+++

# SteepestDescent

It performs steepest descent (gradient descent) method.

It takes one [system]({{<relref "/docs/reference/system">}}) and one [forcefield]({{<relref "/docs/reference/forcefields">}}) to run the simulation.

## Example

```toml
[simulator]
type           = "SteepestDescent"
boundary_type  = "Unlimited"
precision      = "double"
delta          = 1e-4
threshold      = 1e-4
step_limit     = 1_000_000
save_step      = 100
```

## Input Reference

- `type`: String
  - To use MolecularDynamicsSimulator, set `"SteepestDescent"`.
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
  - `"sequencial"`: Simulation runs on single core.
- `delta`: Floating
  - The coefficient of particle movement relative to the force.
- `threshold`: Floating
  - If the maximum displacement of the particles is less than this threshold, it stops the simulation.
- `step_limit`: Integer
  - The limit of total number of steps.
  - It stops if the total number of steps would reach to this limit regardless of the convergence.
- `save_step`: Integer
  - The state of the system will be saved at this interval.
  - The last snapshot will be saved regardless of this number.
