+++
title = "G-JFLangevin"
weight = 3000
+++

# G-JFLangevin

`G-JFLangevin` integrator performs constant temperature simulation according to Langevin equation.

{{<katex display>}}
\begin{aligned}
m\frac{d^2 \bold{r}}{dt^2} &= \bold{f}(\bold{r}) - \alpha\bold{v} + \beta(t) \\
\langle\beta(t)\rangle &= 0 \\
\langle\beta(t)\beta(t')\rangle &= 2\alpha k_B T \delta(t - t')
\end{aligned}
{{</katex>}}


This method is developed in the following paper.

- [Niels Grønbech-Jensen & Oded Farago, (2013) Mol.Phys. 111:8, 983-991](https://doi.org/10.1080/00268976.2012.760055)

## Example

```toml
[simulator]
# ...
integrator.type = "G-JFLangevin"
integrator.alphas = [
    {index = 0, alpha = 1.0},
    {index = 1, alpha = 1.0},
    # ...
]
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - Name of the integrator. Here, it is `"G-JFLangevin"`.
- `alphas`: Array of Tables
  - {{<katex>}}\alpha{{</katex>}} of each particle.
- `remove`: Table (Optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
