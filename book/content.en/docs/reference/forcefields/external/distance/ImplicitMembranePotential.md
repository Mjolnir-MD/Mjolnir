+++
title = "ImplicitMembrane"
weight = 300000
+++

# ImplicitMembranePotential

It is a implicit membrane model that stabilizes the region around the plane.

{{<katex display>}}
U(\mathbf{r}) = \sum_i^N k h_i \tanh\left(\mathrm{bend} * \left(|z_i - z_0| - \frac{\mathrm{thickness}}{2}\right)\right)
{{</katex>}}

## Example

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "ImplicitMembrane"
shape.name     = "AxisAlignedPlane"
shape.axis     = "X"
shape.position = 0.0
shape.margin   = 0.5

bend           = 1.0
thickness      = 4.0
interaction_magnitude = 10.0
parameters  = [
    {index = 0, hydrophobicity = 1.0},
    # ...
]
```

## Input Reference

- `bend`: Floating
  - It modelates the slope between stabilized region and the planner region.
- `thickness`: Floating
  - It determines the region where the stabilized.
- `interaction_magnitude`: Floating
  - It determines the strength of the potential. {{<katex>}} k {{</katex>}} in the above equation.
- `parameters`: Array of Tables
  - `index`: Integer
    - The index of the particle. The index is 0-based.
  - `hydrophobicity`: Floating
    - The hydrophobicity of the particle. {{<katex>}} h {{</katex>}} in the above equation.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
