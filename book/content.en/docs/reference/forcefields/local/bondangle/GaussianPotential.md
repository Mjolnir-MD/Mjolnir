+++
title = "Gaussian"
weight = 200000
+++

# GaussianPotential

The well-known, widely-used gaussian function.

{{<katex display>}}
U(v) = k\exp\left(\frac{-(v - v_0)^2}{2\sigma^2}\right)
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Gaussian"
topology    = "none"
parameters  = [
    {indices = [0, 1, 2], v0 = 1.0, k = -100.0, sigma = 5.0},
    # ...
]
```

## Input Reference

- `k`: Floating
  - It determines the strength of the potential.
  - If positive, the gaussian becomes repulsive.
- `sigma`: Floating
  - It determines the width of the gaussian.
- `v0`: Floating
  - It determines the peak of the gaussian.
- `indices`: Array of Integers (length = 3)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
