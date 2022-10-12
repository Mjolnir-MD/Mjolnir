+++
title = "HarmonicGroove"
weight = 100000
+++

# HarmonicGroovePotential

It is a Harmonic potential between particle and a spatial structure.

{{<<katex dispaly>>}}
U(r) = k(v - v_0)^2
{{</katex>}}

## Example

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "HarmonicGroove"
shape.name     = "AxisAlignedPlane"
shape.axis     = "X"
shape.position = 0.0
shape.margin   = 1.0

parameters = [
    {index = 0, k = 1.0, v0 = 1.0},
    # ...
]
```

## Input Reference

- `k`: Floating
    - It determines the strength of the potential.
- `v0`: Floating
    - It determines the native length.
