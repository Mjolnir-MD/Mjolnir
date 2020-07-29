+++
title  = "SimulatedAnnealing"
weight = 2000
+++

# SimulatedAnnealing

It runs simulated annealing simulation.

It takes one [system]({{<relref "/docs/reference/system">}}) and one [forcefield]({{<relref "/docs/reference/forcefields">}}) to run the simulation.

## Example

```toml
[simulator]
type           = "SimulatedAnnealing"
boundary_type  = "Unlimited"
precision      = "double"
parallelism    = "OpenMP" # optional
delta_t        = 0.1
total_step     = 50_000
save_step      = 100
each_step      = 100
schedule.type  = "linear"
schedule.begin = 300.0 # temperature in [K]
schedule.end   = 150.0 # temperature in [K]

integrator.type = "UnderdampedLangevin"
integrator.parameters = [
    # ...
]
```

## Input Reference

- `type`: String
  - Name of the simulator. Here, it is `"SimulatedAnnealing"`.
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
- `each_step`: Integer
  - The temperature of the system will be updated at this interval.
- `schedule`: Table
  - The schedule for temperature control. See below.
- `integrator`: Table
  - The time integration method to be used.
  - ["BAOABLangevin"]({{<relref "/docs/reference/integrators/BAOABLangevinIntegrator.md">}})
  - ["g-BAOABLangevin"]({{<relref "/docs/reference/integrators/gBAOABLangevinIntegrator.md">}})
  - ["UnderdampedLangevin"]({{<relref "/docs/reference/integrators/UnderdampedLangevinIntegrator.md">}})
  - For detail, see [`integrators`]({{<relref "/docs/reference/integrators">}}).

### `schedule` Table

- `type`: String
  - Type of the control curve. Available curves are below.
  - `"linear"`
- `begin`: Floating
  - The temperature at the first frame. The unit is [K].
- `end`: Floating
  - The temperature at the last frame. The unit is [K].
