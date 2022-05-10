+++
title = "3SPN2Bond"
weight = 1000000
+++

# 3SPN2BondPotential

It is a component of 3SPN2 potential.

{{<katex display>}}
U(v) = k(v - v_0)^2 + 100k (v - v_0)^4
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "3SPN2Bond"
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100, v0 = 1.0, k = 10.0},
    # ...
]
```

## Input Reference

- `k`: Floating
  - It determines the strength of the potential.
- `v0`: Floating
  - It determines the native length.
- `indices`: Array of Integers (length = 2)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
