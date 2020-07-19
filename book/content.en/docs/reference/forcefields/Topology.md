+++
title = "Topology"
weight = 6000
+++

# Topology

In some cases, `LocalForceField` affects `GlobalForceField`.
For example, excluded volume interaction would not be applied to a pair of particles that interact with each other via bond length interaction.

To share the information among forcefields, `Topology` is constructed.

It contains graph structure where particles are represented as nodes and local interactions are represented as edges.

[`LocalForceField`]({{<relref local>}}) can name the interaction.
[`GlobalForceField`]({{<relref global>}}) can ignore specific local interactions.

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    {indices = [0, 1], ... },
    # ...
]

[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"
ignore.particles_within.bond    = 3 # ignore particles within 3 bonds.
ignore.particles_within.contact = 1 # ignore particles within 1 contact.
# ...
parameters = [
    # ...
]
```
