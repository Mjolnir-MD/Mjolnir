+++
title = "GoContact"
weight = 300000
+++

# GoContactPotential

10-12 Lennard Jones potential.

{{<katex display>}}
U(r) = k\left[5\left(\frac{r_0}{r}\right)^{12} - 6\left(\frac{r_0}{r}\right)^{10}\right]
{{</katex>}}

It is mainly used in a structure-based Coarse-Grained models to represent attractive pairs of particles.

## Example

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "GoContact"
topology    = "contact"
parameters = [
    {indices = [0, 1], v0 = 1.0, k = 0.1},
    # ...
]
```

## Input Reference

- `v0`: Floating
  - The native distance. {{<katex>}} r_0 {{</katex>}} in the above equation.
- `k`: Floating
  - The strength of the potential.
- `indices`: Array of Integers (length = 2)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
