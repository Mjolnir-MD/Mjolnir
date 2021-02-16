+++
title = "g-BAOABLangevin"
weight = 2000
+++

# g-BAOABLangevin

`g-BAOABLangevin` integrator performs constant temperature simulation according to Langevin equation.

Unlike `BAOABLangevin`, it handles bond length constraints appropreately.

`g-BAOABLangevin` is developed in the following paper.

- [Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016) 472: 20160138](https://doi.org/10.1098/rspa.2016.0138)

## Example

```toml
[simulator]
integrator.type = "g-BAOABLangevin"
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - Name of the integrator. Here, it is `"g-BAOABLangevin"`.
- `gammas`: Array of Tables
  - {{<katex>}}\gamma_i{{</katex>}} of the particles.
- `remove`: Table (optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
