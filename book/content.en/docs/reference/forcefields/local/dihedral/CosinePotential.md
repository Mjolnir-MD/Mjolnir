+++
title = "Cosine"
weight = 100000
+++

# CosinePotential

Cosine potential is a periodic alternative of harmonic potential for dihedral angle interaction.

The functional form is the following.

{{<katex display>}}
U(v) = k\left(1 + \cos(n(v - v_0))\right)
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Cosine"
topology    = "none"
parameters = [
    {indices = [0, 1, 2, 3], v0 = 1.57, k = 10.0, n = 1},
]
```

## Input Reference

- `k`: Floating
  - It determines the strength of the potential.
- `n`: Integer
  - It determines the number of minima in the potential.
- `v0`: Floating
  - It determines the position of minima.
- `indices`: Array of Integers (length = 4)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
