+++
title = "StoichiometricUniformCubicPan"
weight = 100000
+++

# StoichiometricUniformCubicPanPotential

The simple cubic functional potential for `StoichiometricInteraction` which represent pan shape.

{{<katex display>}}
U(r) = \begin{cases}
- 1 & 0 < r \leq d \\
- \left(1 - 3 \left(\frac{r - d}{\Delta r} \right)^2 + 2 \left( \frac{r - d}{\Delta r} \right)^3 \right) & d < r \leq d + \Delta r\\
0 & d + \Delta r < r
\end{cases}
{{</katex>}}

This is for the case where all particles has the same parameter.

## Example

```toml
[[forcefields.global]]
interaction              = "Stoichiometric"
potential                = "StoichiometricUniformCubicPan"
ignore.molecule          = "Nothing"
spatial_partition.type   = {type = "CellList", margin = 0.2}
spatial_partition.margin = 0.2
particle_kinds = [
    {name = "A", coef = 1},
    {name = "B", coef = 1}
]
epsilon = 5.0  # parameter for StoichiometricInteraction
v0      = 10.0 # parameter for StoichiometricUniformCubicPanPotential
range   = 10.0 # parameter for StoichiometricUniformCubicPanPotential
paramters = [
    {index = 0, kind = "A"},
    # ...
]
```

## Input Reference
- `index`: Integer
  - The index of the particle.
- `offset`: Integer (optional. By default, 0)
  - Offset of the index.
- `v0`: Floating
  - It determines the range of the bottom of pan.
- `range`: Floating
  - It determines the range of the edge of pan.
