# BAOABLangevin

`BAOABLangevin` integrator runs Langevin dynamics.

This method is developed in the following paper.

- Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)

## Example

```toml
[simulator]
integrator.type = "BAOABLangevin"
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[Simulator]`](Simulator.md) table.

- `type`: String
  - The name of [Integrator](Integrator.md). Here, it is `"BAOABLangevin"`.
- `remove`: Table (optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
- `gammas`: Array of Tables
  - $$\gamma_i$$ of the particles.
