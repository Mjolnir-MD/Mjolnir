+++
title = "GFWNPTLangevin"
weight = 4000
+++

# GFWNPTLangevin

`GFWNPTLangevin` integrator performs constant temperature and constant pressure simulation according to Langevin equation.

This method is introduced in the following paper.

- [Xingyu Gao, Jun Fang, and Han Wang. J. Chem. Phys. (2016) 144, 124113](https://doi.org/10.1063/1.4944909)

## Example

```toml
[simulator]
# ...
integrator.type = "GFWNPTLangevin"
integrator.chi  = 0.0
integrator.cell_mass  = [1e3, 1e3, 1e3]
integrator.cell_gamma = [0.1, 0.1, 0.1]
integrator.gammas = [
    {index = 0, gamma = 0.1},
    {index = 1, gamma = 0.1},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - Name of the integrator. Here, it is `"GFWNPTLangevin"`.
- `chi`: Floating
  - {{<katex>}} chi {{</katex>}} value introduced in the paper.
- `cell_mass`: Array of Floating numbers
  - virtual mass of the box in each direction.
- `cell_gamma`: Array of Floating numbers
  - friction coefficient of the cell.
- `gammas`: Array of Tables
  - {{<katex>}}\gamma_i{{</katex>}} of the particles.
- `remove`: Table (Optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
