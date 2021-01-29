+++
title = "UnderdampedLangevin"
weight = 4000
+++

# UnderdampedLangevin

`UnderdampedLangevin` integrator performs constant temperature simulation according to Langevin equation.

This method is developed in the following paper.

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

It is the same method that is employed in CafeMol.

## Example

```toml
[simulator]
# ...
integrator.type = "UnderdampedLangevin"
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[Simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - The name of Integrator. Here, it is `"UnderdampedLangevin"`.
- `gammas`: Array of Tables
  - {{<katex>}}\gamma_i{{</katex>}} of the particles.
- `remove`: Table (Optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
