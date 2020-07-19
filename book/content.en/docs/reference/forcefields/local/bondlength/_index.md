+++
title = "BondLength"
weight = 10000
bookCollapseSection = true
+++

# BondLengthInteraction

`BondLengthInteraction` is an interaction that depends on the distance between particles.

The following potentials are available.
- [Harmonic]({{<relref "HarmonicPotential.md">}})
- [Gaussian]({{<relref "GaussianPotential.md">}})
- [GoContact]({{<relref "GoContactPotential.md">}})
- [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
- [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})
- [WormLikeChain]({{<relref "WormLikeChainPotential.md">}})
- [3SPN2Bond]({{<relref "3SPN2BondPotential.md">}})

## Example

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1], offset = 100, ... },
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"BondLength"`.
- `potential`: String
  - The following potentials are available.
  - [Harmonic]({{<relref "HarmonicPotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [GoContact]({{<relref "GoContactPotential.md">}})
  - [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
  - [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})
  - [WormLikeChain]({{<relref "WormLikeChainPotential.md">}})
  - [3SPN2Bond]({{<relref "3SPN2BondPotential.md">}})
- `topology`: String
  - Name of the connection in [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 2)
    - Indices of particles that interact with each other. The index is 0-based.
  - `offset`: Integer(Optional. By default, 0.)
    - Offset of index.
  - The other parameter depends on the specified potential.
