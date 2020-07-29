+++
title = "Pair"
weight = 10000
bookCollapseSection = true
+++

# GlobalPair

`GlobalPairInteraction` will be applied to all the possible pairs of particles participating in the interaction.

The following potentials are available.

- [LennardJones]({{<relref "LennardJonesPotential.md">}})
- [UniformLennardJones]({{<relref "UniformLennardJonesPotential.md">}})
- [DebyeHuckel](DebyeHuckelPotential.md)
- [ExcludedVolume](ExcludedVolumePotential.md)
- [InversePower](InversePowerPotential.md)
- [HardCoreExcludedVolume](HardCoreExcludedVolumePotential.md)
- [3SPN2ExcludedVolume](3SPN2ExcludedVolumePotential.md)

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
ignore.groups.intra = ["chain-A"]
ignore.groups.inter = [["chain-B", "chain-C"]]
spatial_partition.type = "CellList"
spatial_partition.margin = 0.2
parameters = [
    {index = 0, offset = 100, ...}, # required parameter depends on the potential.
    # ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"Pair"`.
- `potential`: String
  - The following potentials are available.
  - [`"LennardJones"`]({{<relref "LennardJonesPotential.md">}})
  - [`"UniformLennardJones"`]({{<relref "UniformLennardJonesPotential.md">}})
  - [`"DebyeHuckel"`](DebyeHuckelPotential.md)
  - [`"ExcludedVolume"`](ExcludedVolumePotential.md)
  - [`"InversePower"`](InversePowerPotential.md)
  - [`"HardCoreExcludedVolume"`](HardCoreExcludedVolumePotential.md)
- `ignore`: Table
  - It describes the condition when the pair of particles does not interact to each other.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#ignore">}})
- `spatial_partition`: Table
  - It specifies the algorithm to construct a neighbor list.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#spatial_partition">}})
- `parameters`: Array of Tables
  - `index`: Integer
    - The index of a particle. The index is 0-based.
  - `offset`: Integer (Optional. By default, 0.)
    - The offset value for the index.
