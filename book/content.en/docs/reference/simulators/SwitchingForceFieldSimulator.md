+++
title = "SwitchingForceField"
weight = 4000
+++

# SwitchingForceField

It performs molecular dynamics simulation and changes forcefields at the given frame.

It takes one [system]({{<relref "/docs/reference/system">}}) and several [forcefields]({{<relref "/docs/reference/forcefields">}}) to run the simulation.

## Example

```toml
[simulator]
type          = "SwitchingForceField"
boundary_type = "Unlimited"
precision     = "double"
delta_t       = 0.1
total_step    = 3_000_000
save_step     = 100
seed          = 2374

integrator.type = "BAOABLangevin"
integrator.parameters = [
{index =  0, gamma = 1.00},
{index =  1, gamma = 1.00},
]

schedule      = [
    {until = 1_000_000, forcefield = "close"},
    {until = 2_000_000, forcefield = "open"},
    {until = 3_000_000, forcefield = "close"},
]

[[forcefields]]
name = "close"
[[forcefields.local]]
# ...

[[forcefields]]
name = "open"
[[forcefields.local]]
# ...
```

## Input Reference

- `type`: String
  - Name of the simulator. Here, it is `"SwitchingForceField"`.
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
- `schedule`: Array of Tables
  - The schedule for forcefield control. See below.

### `schedule`

`schedule` is an array of tables and each table has the following fields.

- `until`: Integer
  - It uses the specified forcefield until this frame.
- `forcefield`: String
  - Name of the forcefield.

```toml
schedule      = [
    {until = 1_000_000, forcefield = "close"},
    {until = 2_000_000, forcefield = "open"},
    {until = 3_000_000, forcefield = "close"},
]
```

The name of the `forcefield` is defined via `name` under `[[forcefields]]`.

```toml
[[forcefields]]
name = "close"

[[forcefields.local]]
# This is a component of "close" forcefield.
interaction = "BondLength"
potential = "Harmonic"
# ...

[[forcefields.local]]
# This is a component of "close" forcefield.
interaction = "BondAngle"
potential = "Harmonic"
# ...

[[forcefields]]
name = "open"

[[forcefields.local]]
# This is a component of "open" forcefield.
# ...
```
