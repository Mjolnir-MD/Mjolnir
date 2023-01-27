+++
title = "FlexibleLocal"
weight = 300000
+++

# FlexibleLocalDihedral

It represents dihedral angle distribution of flexible regions of proteins in a Coarse-Grained protein model.

{{<katex display>}}
U(\phi) = C + \sum_{n=1}^{3}\left( k_n^{\sin} \sin(n\phi) + k_n^{\cos} \cos(n\phi)\right)
{{</katex>}}

It is a component of AICG2+ Coarse-Grained protein model.

It is developed in the following paper.

- T. Terakawa and S. Takada, (2011) Biophys J 

## Example

```toml
[[forcefields.local]]
parameters = [
    {indices = [0,1,2,3], k = 1.0, coef = [2.2356,  0.4119, -0.1283,  0.0229, -0.2708, -0.0085, -0.0641]},
    # ...
]
```

## Input Reference

- `k`: Floating
  - It determines the strength of the potential.
- `coef`: Array of Floatings (length = 7)
  - {{<katex>}} C, k_1^{\cos}, k_1^{\sin}, k_2^{\cos}, k_2^{\sin}, k_3^{\cos}, k_3^{\sin} {{</katex>}}, respectively.
- `indices`: Array of Integers (length = 4)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
