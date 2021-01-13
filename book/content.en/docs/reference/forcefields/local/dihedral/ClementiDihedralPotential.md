+++
title = "ClementiDihedral"
weight = 400000
+++

# ClementiDihedral

It is specific to Off-lattice Go forcefield.

{{<katex display>}}
U(v) = k_1(1-\cos(v-v_0)) + k_3(1-\cos(3(v-v_0)))
{{</katex>}}

It is introduced in the following paper.

- C. Clementi, H. Nymeyer, J. Onuchic, (2000) JMB

## Example

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "ClementiDihedral"
topology    = "none"
parameters  = [
    {indices = [0,1,2,3], v0 = -2.20, k1 = 1.0, k3 = 0.5},
    # ...
]
```

## Input Reference

- `v0`: Floating
  - It determines the minimum of the potenital.
- `k1`: Floating
  - It determines the strength of the potential.
- `k3`: Floating
  - It determines the strength of the potential.
- `indices`: Array of Integers (length = 4)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
