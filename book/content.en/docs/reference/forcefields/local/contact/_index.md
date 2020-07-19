+++
title = "Contact"
weight = 40000
bookCollapseSection = true
+++

# ContactInteraction

`ContactInteraction` is an interaction that depends on the distance between particles.

Basically it is the same as `BondLength`.
But it contains a neighboring list and it skips the calculation if the distance exceeds the cutoff.

The following potentials are available.

- [Gaussian]({{<relref "GaussianPotential.md">}})
- [GoContact]({{<relref "GoContactPotential.md">}})
- [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
- [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})


## Example

```toml
[[forcefields.local]]
interaction = "Contact"
potential   = "GoContact"
topology    = "bond"
margin      = 0.5 # relative length to longest cutoff
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1], ... },
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"BondLength"`.
- `potential`: String
  - The following potentials are available.
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [GoContact]({{<relref "GoContactPotential.md">}})
  - [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
  - [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})
- `topology`: String
  - Name of the connection in [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).
- `margin`: Floating
  - The margin in the neighboring list.
  - It is relative to the maximum cutoff distance.
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 2)
    - Indices of particles that interact with each other. The index is 0-based.
  - `offset`: Integer(Optional. By default, 0.)
    - Offset of index.
  - The other parameter depends on the specified potential.
