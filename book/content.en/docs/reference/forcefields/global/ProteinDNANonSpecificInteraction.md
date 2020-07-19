+++
title = "PDNS"
weight = 30000
+++

# ProteinDNANonSpecificInteraction

`ProteinDNANonSpecificInteraction` is a coarse-grained model of hydrogen bond, especially formed between protein and DNA.

It is developed in the following paper.

- T.Niina‡, G.B.Brandani‡, C.Tan and S.Takada, (2017) PLoS Comput Biol. (‡: co-1st)

## Example

```toml
[[forcefields.global]]
interaction = "PDNS"
potential   = "PDNS"
spatial_partition.type   = "VerletList"
spatial_partition.margin = 0.5
sigma = 1.0
delta = 0.17453
parameters  = [
{index =    2, kind = "DNA", S3 = 1},
{index =    5, kind = "DNA", S3 = 4},
# ...
{index = 1000, offset = 100, kind = "Protein", PN =  999, PC = 1001, k = -1.2, r0 = 5.0, theta0 = 1.57, phi0 = 1.73},
{index = 1023, offset = 100, kind = "Protein", PN = 1022, PC = 1024, k = -1.2, r0 = 6.0, theta0 = 1.57, phi0 = 1.73},
# ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"PDNS"`.
- `potential`: String
  - Name of the potential. Here, it is `"PDNS"`.
- `sigma`: Floating
  - The width of the attractive potential along the radial direction.
- `delta`: Floating
  - The width of the attractive potential along the angle direction.
- `parameters`: Array of Tables
  - `index`: Integer
    - The index of the particle.
    - In `DNA` case, the index of the phosphate.
  - `offset`: Integer (Optional. By default, 0)
    - Offset of the index.
  - `kind`: String
    - There are 2 kinds of particles, `DNA` and `Protein`.
  - `S3`: Integer
    - Only required for `DNA`. The index of sugar particle of adjacent nucleotide towords the 3' end.
  - `PN, PC`: Integer
    - Only required for `Protein`. The adjacent particle towords N- and C-terminus, respectively.
  - `k`: Floating
    - The strength of the potential.
  - `r0`, `theta0`, `phi0`: Floating
    - The native conformation.

