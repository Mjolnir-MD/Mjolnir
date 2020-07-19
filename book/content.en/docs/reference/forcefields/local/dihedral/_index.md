+++
title = "DihedralAngle"
weight = 30000
bookCollapseSection = true
+++

# DihedralAngleInteraction

`DihedralAngleInteraction` depends on the angle formed by 2 planes, where each plane is defined by 3 particles.

The following potentials are available.

- [Cosine]({{<relref "CosinePotential.md">}})
- [Gaussian]({{<relref "GaussianPotential.md">}})
- [ClementiDihedral]({{<relref "ClementiDihedralPotential.md">}})
- [FlexibleLocalDihedral]({{<relref "FlexibleLocalDihedralPotential.md">}})

## Example

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Cosine"
topology    = "none"
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1, 2, 3], ... },
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"BondLength"`.
- `potential`: String
  - The following potentials are available.
  - [Cosine]({{<relref "CosinePotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [ClementiDihedral]({{<relref "ClementiDihedralPotential.md">}})
  - [FlexibleLocalDihedral]({{<relref "FlexibleLocalDihedralPotential.md">}})
- `topology`: String
  - Name of the connection in [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 3)
    - Indices of particles that interact with each other. The index is 0-based.
  - `offset`: Integer(Optional. By default, 0.)
    - Offset of index.
  - The other parameter depends on the specified potential.

## Combination

Some forcefields apply different potential functions to the same set of particles.
Since dihedral angle calculation is costly operation relative to other local interactions, it is better to use those functions at once.
The following frequently-used combination is supported.

- `"Gaussian+FlexibleLocalDihedral"`
- `"Gaussian+Cosine"`

Parameters of each potential is defined as an inline table, like: `PotentialName = {}`.
For example,

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Gaussian+FlexibleLocalDihedral"
topology    = "none"
env.ALA-ALA = [2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375]
parameters  = [
{indices = [0,1,2,3], Gaussian = {v0=-2.2, k=-0.43,sigma=0.15}, FlexibleLocalDihedral = {k=1.0, coef="ALA-ALA"}},
# ...
]
```
