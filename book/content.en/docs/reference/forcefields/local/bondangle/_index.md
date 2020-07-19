+++
title = "BondAngle"
weight = 20000
bookCollapseSection = true
+++

# BondAngleInteraction

`BondAngleInteraction` depends on the angle formed by 3 particles.

The following potentials are available.
- [Harmonic]({{<relref "HarmonicPotential.md">}})
- [Gaussian]({{<relref "GaussianPotential.md">}})
- [FlexibleLocalAngle]({{<relref "FlexibleLocalAnglePotential.md">}})

## Example

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Harmonic"
topology    = "none"
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1, 2], offset = 100, ... },
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"BondLength"`.
- `potential`: String
  - The following potentials are available.
  - [Harmonic]({{<relref "HarmonicPotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [FlexibleLocalAngle]({{<relref "FlexibleLocalAnglePotential.md">}})
- `topology`: String
  - Name of the connection in [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 3)
    - Indices of particles that interact with each other. The index is 0-based.
  - `offset`: Integer(Optional. By default, 0.)
    - Offset of index.
  - The other parameter depends on the specified potential.
