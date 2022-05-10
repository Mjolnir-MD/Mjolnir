+++
title = "RepulsiveMBasin"
weight = 600000
+++

# RepulsiveMBasinContactPotential

Repulsive part of the [Go-Contact potential]({{<relref "GoContactPotential.md">}}) specifically used in MultipleBasin forcefield.

{{<katex display>}}
U(r) = k\max\left(0, 5\left(\frac{r_0}{r}\right)^{12} - 6\left(\frac{r_0}{r}\right)^{10})
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "Contact"
potential   = "AttractiveMBasinContact"
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
