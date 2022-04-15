+++
title = "UniformCubicPan"
weight = 200000
+++

# UniformCubicPanPotential

The simple cubic functional potential which represent pan shape.

{{<katex display>}}
U(r) = \begin{cases}
-\varepsilon & 0 < r \leq d \\
-\varepsilon \left(1 - 3 \left(\frac{r - d}{\Delta r} \right)^2 + 2 \left( \frac{r - d}{\Delta r} \right)^3 \right) & d < r \leq d + \Delta r\\
0 & d + \Delta r < r
\end{cases}
{{</katex>}}

This is for the case where all particles has the same parameter.

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "CubicPanPotential"
ignore.molecule = "Nothing"
ignore.partcles_within.bond     = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

epsilon = 1.0
v0      = 9.0
range   = 10.0
parameters = [
    {index = 0, offset = 100}, # to control which particle participantes
]
```

## Input Reference

- `index`: Integer
  - The index of the particle.
- `offset`: Integer (optional. By default, 0)
  - Offset of the index.
- `epsilon`: Floating
  - It determines the strength of the potential.
- `v0`: Floating
  - It determines the range of the bottom of pan.
- `range`: Floating
  - It determines the range of the edge of pan.
