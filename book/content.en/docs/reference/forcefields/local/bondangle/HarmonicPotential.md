+++
title = "Harmonic"
weight = 100000
+++

# HarmonicPotential

Well-known and widely-used harmonic function.

{{<katex display>}}
U(v) = k(v-v_0)^2
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    {indices = [0, 1, 2], offset = 100, v0 = 1.0, k = 100.0},
    # ...
]
```

## Input Reference

- `k`: Floating
  - It determines the strength of the potential.
- `v0`: Floating
  - It determines the native length.
- `indices`: Array of Integers (length = 3)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
