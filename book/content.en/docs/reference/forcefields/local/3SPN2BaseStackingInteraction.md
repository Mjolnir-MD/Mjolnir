+++
title = "3SPN2 BaseStacking"
weight = 50000
+++

# 3SPN2BaseStackingInteraction

`3SPN2BaseStackingInteraction` is specific to 3SPN2 Coarse-Grained DNA model.

The following potentials are available.
- `"3SPN2"` (Hinckley et al., (2013) JCP)
- `"3SPN2C"` (Freeman et al., (2014) JCP)

## Example

```toml
[[forcefields.local]]
interaction = "3SPN2BaseStacking"
potential   = "3SPN2"
topology    = "nucleotide"
parameters  = [
# `nucleotide` index starts from 5' and ends at 3'.
{strand = 0, nucleotide =  0,          S =   0, B =   1, offset = 100, Base = "A"},
{strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, offset = 100, Base = "T"},
{strand = 0, nucleotide =  2, P =   5, S =   6, B =   7, offset = 100, Base = "C"},
{strand = 0, nucleotide =  3, P =   8, S =   9, B =  10, offset = 100, Base = "G"},
{strand = 1, nucleotide =  4,          S =  11, B =  12, offset = 100, Base = "C"},
{strand = 1, nucleotide =  5, P =  13, S =  14, B =  15, offset = 100, Base = "G"},
{strand = 1, nucleotide =  6, P =  16, S =  17, B =  18, offset = 100, Base = "A"},
{strand = 1, nucleotide =  7, P =  19, S =  20, B =  21, offset = 100, Base = "T"},
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"3SPN2BaseStacking"`.
- `potential`: String
  - Name of the potential. One of the following.
  - `3SPN2`
  - `3SPN2C`
- `topology`: String
  - Name of the topology. Debye-Huckel and 3SPN2BaseBaseInteraction should handle this correctly.
  - This interaction adds edges between all the pair of nucleotide particles that are next to each other.
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
