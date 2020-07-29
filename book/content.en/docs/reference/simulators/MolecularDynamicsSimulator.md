+++
title  = "MolecularDynamics"
weight = 1000
+++

# MolecularDynamics

It performs a normal molecular dynamics simulation.

It takes one [system]({{<relref "/docs/reference/system">}}) and one [forcefield]({{<relref "/docs/reference/forcefields">}}) to run the simulation.

## Example

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "PeriodicCuboid"
precision     = "double"
parallelism   = "OpenMP" # optional
seed          = 12345
delta_t       = 0.1
total_step    = 50_000
save_step     = 100
integrator.type = "VelocityVerlet"
```

## Input Reference

- `type`: String
  - To use MolecularDynamicsSimulator, set `"MolecularDynamics"`.
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
- `seed`: Integer
  - Random number generator will be initialized by this value.
- `delta_t`: Floating
  - Time step of the simulation. The unit depends on the unit system defined in [`[units]`]({{<relref "/docs/reference/units">}})
- `total_step`: Integer
  - Total time step of the simulation.
- `save_step`: Integer
  - The state of the system will be saved at this interval.
- `integrator`: Table
  - The time integration method to be used.
  - ["BAOABLangevin"]({{<relref "/docs/reference/integrators/BAOABLangevinIntegrator.md">}})
  - ["g-BAOABLangevin"]({{<relref "/docs/reference/integrators/gBAOABLangevinIntegrator.md">}})
  - ["UnderdampedLangevin"]({{<relref "/docs/reference/integrators/UnderdampedLangevinIntegrator.md">}})
  - ["VelocityVerlet"]({{<relref "/docs/reference/integrators/VelocityVerletIntegrator.md">}})
  - For detail, see [`integrators`]({{<relref "/docs/reference/integrators">}}).
