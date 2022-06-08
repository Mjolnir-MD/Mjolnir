+++
title = "Stoichiometric"
weight = 20000
+++

# StoichiometricInteraction

`StoichiometricInteraction` is the interaction for representing the strict valency of binding between particles. This interaction applies heterotypic pairs defined in `parameter_kinds`. You can specify stoichiometric coefficient of the interaction.

{{<katex display>}}
V(r_{ij}) = - \epsilon \sum_{i \in a, j \in b} \frac{U(r_{ij})}{\max({\rm CO}_i, {\rm CO}_j, 1)}
{{</katex>}}

where {{<katex>}}{\rm CO}_i{{</katex>}} and {{<katex>}}{\rm CO}_j{{</katex>}} are

{{<katex display>}}
{\rm CO}_i = \frac{\sum_{j \in b} u(r_{ij})}{m}, {\rm CO}_j = \frac{\sum_{i \in a} u(r_{ij})}{n}
{{</katex>}}

.  {{<katex>}}m{{</katex>}} means the stoichiometric coefficient of b particle and {{<katex>}}n{{</katex>}} means the stoichiometric coefficient of a particle.


## Example

```toml
[[forcefields.global]]
interaction              = "Stoichiometric"
potential                = "StoichiometricUniformCubicPan"
ignore.molecule          = "Nothing"
spatial_partition.type   = {type = "CellLIst", margin = 0.2}
sparial_partition.margin = 2.9100810346402386
parameter_kinds = [
    {name = "A", coef = 1},
    {name = "B", coef = 1}
]
v0    = 10.0 # parameter for StoichiometricUniformCubicPanPotential
range = 10.0 # parameter for StoichiometricUniformCubicPanPotential
parameters = [
{index = 0, name = "A"},
{index = 1, name = "A"},
{index = 2, name = "B"},
{index = 3, name = "B"}
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is "Stoichiometric"
- `potential` : String
  - The fllowing potential is available.
    - [`"StoichiometricUniformCubicPan"`]({{<relref "StoichiometricUniformCubicPanPotential.md">}})
- `ignore`: Table
  - It describes the condition when the pair of particles does not interact to each other.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#ignore">}}).
- `spatial_partition`: Table
  - It specifies the algorithm to construct a neighbor list.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#ignore">}}).
- `particle_kinds`: Array of Tables
  - It describes the names and stoichiometric coefficients of particles in this interaction. This interaction only allow two kind case, so *the length of this array must be two*.
  - `name`: String
    - It determin the particle name for the parameters section.
  - `coef`: Integer
    - It determin the stoichiometric coefficients of the particle kinds.
- `parameters`: Array of Tables
  - `index`: Integer
    - Index of the particle.
  - `offset`: Integer (Optional. By default, 0.)
    - The offset value for the index.
  - `name`: String
    - The name of particle determined in `particle_kinds` section.
