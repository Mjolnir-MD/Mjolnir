+++
title  = "Simulators"
weight = 300
bookCollapseSection = true
+++

# `[simulator]`

`[simulator]` table defines what kind of simulation it will perform.

## Input Reference

### Common Part

Some of the fields are common to all the simulators.

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "Periodic"
precision     = "double"
parallelism   = "OpenMP" # optional
```

- `type`: String
  - It defines simulation method. E.g. simulated annealing, steepest descent, etc.
  - See the following section, "Available Simulators".
- `boundary_type`: String
  - `"Unlimited"`: No boundary will be applied.
  - `"Periodic"`: It applies periodic boundary of which shape is recutangular box.
  - `"PeriodicCuboid"`: The same as the previous option.
- `precition`: String
  - `"float"`: Use 32-bit floating point type.
  - `"double"`: Use 64-bit floating point type.
- `parallelism`: String (Optional. By default, `"sequencial"`.)
  - `"OpenMP"`: Use OpenMP implementation.
  - `"sequencial"`: Run on single core.
- `forcefield`: Table (Optional. By default, none.)
  - For detail, see [MultipleBasinForceField]({{<relref "/docs/reference/forcefields/MultipleBasinForceField.md">}}).

### Available Simulators

- [MolecularDynamics]({{<relref "MolecularDynamicsSimulator.md">}})
  - It performs normal molecular dynamics simulation.
- [SimulatedAnnealing]({{<relref "SimulatedAnnealingSimulator.md">}})
  - It performs [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing) simulation with given forcefield.
- [SteepestDescent]({{<relref "SteepestDescentSimulator.md">}})
  - It performs [steepest descent method](https://en.wikipedia.org/wiki/Gradient_descent) with given forcefield.
- [SwitchingForceField]({{<relref "SwitchingForceFieldSimulator.md">}})
  - It performs normal molecular dynamics simulation but changes forcefields as scheduled order.
- [EnergyCalculation]({{<relref "EnergyCalculationSimulator.md">}})
  - It does not perform any simulation but calculates energy from trajectory with given forcefield.
