+++
title = "Stoichiometric"
weight = 20000
+++

# StoichiometricInteraction

`StoichiometricInteraction` is the interaction for represent the strict valency of binding between particles. You can specify stoichiometric coefficient  of the interaction by `parameter_kinds` table.


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
    - `StoichiometricUniformCubicPan`
- `ignore`: Table
  - It describes the condition when the pair of particles does not interact to each other.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#ignore">}}).
- `spatial_partition`: Table
  - It specifies the algorithm to construct a neighbor list.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#ignore">}}).
- `parameter_kinds`: Array of Tables
  - It describes the names and stoichiometric coefficients of particles in this interaction. This interaction only allow two kind case, so *the length of this array must be two*.
  - `name`: String
    - It determin the particle name for the parameters section.
  - `coef`: Integer
    - It determin the stoichiometric coefficients of the particle kinds.
- `parameters`: Array of Tables
  - `strand`: Integer
    - Index of the particle.
  - `name`: String
    - The name of particle determined in `parameter_kinds` section.
