+++
title = "3SPN2Base"
weight = 20000
+++

# 3SPN2BaseBaseInteraction

`3SPN2BaseBaseInteraction` is speicfic to the 3SPN2 Coarse-Grained DNA model.

The following potentials are available.

- `"3SPN2"`: Hinckley et al., (2013) JCP
- `"3SPN2C"`: Freeman et al., (2014) JCP

## Example

```toml
[[forcefields.global]]
interaction = "3SPN2BaseBase"
potential   = "3SPN2"
ignore.particles_within.nucleotide = 3
spatial_partition = {type = "CellList", margin = 0.2}
parameters  = [
# `nucleotide` index starts from 5' and ends at 3'.
{strand = 0, nucleotide =  0,          S =   0, B =   1, offset = 100, Base = "A"},
{strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, offset = 100, Base = "T"},
{strand = 0, nucleotide =  2, P =   5, S =   6, B =   7, offset = 100, Base = "C"},
# ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"3SPN2BaseBase"`.
- `potential`: String
  - The following potentials are available.
  - `3SPN2`
  - `3SPN2C`
- `ignore`: Table
  - It describes the condition when the pair of particles does not interact to each other.
  - Since 3SPN2 only allows base pairs between nucleotides that are at least 3 nucleotides distant. Set `nucleotide = 3` and name [`3SPN2BaseStacking.topology`]({{<relref "/docs/reference/forcefields/local/3SPN2BaseStackingInteraction.md">}}) as `nucleotide`.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#ignore">}})
- `spatial_partition`: Table
  - It specifies the algorithm to construct a neighbor list.
  - For detail, see [the ignore section of GlobalForceField]({{<relref "/docs/reference/forcefields/global#spatial_partition">}})
- `parameters`: Array of Tables
  - `strand`: Integer
    - Index of the strand.
  - `nucleotide`: Integer
    - Index of the nucleotide.
  - `P, S, B`: Integer
    - Indices of particles that correspond to phosphate (P), sugar (S), and base (B).
    - The index is 0-based.
  - `offset`: Integer (Optional. By default, 0.)
    - The offset value for the index.
  - `Base`: String
    - One of `"A"`, `"T"`, `"C"` or `"G"`.
  - Normally nucleotide at the edge of the DNA does not have phosphate.
