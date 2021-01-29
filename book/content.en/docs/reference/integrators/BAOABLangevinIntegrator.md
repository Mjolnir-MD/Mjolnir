+++
title = "BAOABLangevin"
weight = 1000
+++

# BAOABLangevin

`BAOABLangevin` integrator performs constant temperature simulation according to Langevin equation.

This method is developed in the following paper.

- [Benedict Leimkuhler and Charles Matthews. Appl. Math. Res. Exp. (2013) 2013:1, pp. 34-56](https://doi.org/10.1093/amrx/abs010)
- [Benedict Leimkuhler and Charles Matthews. J. Chem. Phys. (2013) 138:17, 174102](https://doi.org/10.1063/1.4802990)

## Example

```toml
[simulator]
# ...
integrator.type = "BAOABLangevin"
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - Name of the integrator. Here, it is `"BAOABLangevin"`.
- `gammas`: Array of Tables
  - {{<katex>}}\gamma_i{{</katex>}} of the particles.
- `remove`: Table (Optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
