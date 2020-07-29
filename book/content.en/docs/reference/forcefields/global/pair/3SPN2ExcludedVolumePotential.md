+++
title = "3SPN2ExcludedVolume"
weight = 700000
+++

# 3SPN2ExcludedVolumePotential

It is specific to 3SPN2 coarse-grained DNA model.

{{<katex display>}}
U(r) =
\begin{cases}
\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - 2\left(\frac{\sigma_{ij}}{r}\right)^6\right] + \epsilon & (r < \sigma_{ij})\\
0 & (r \geq \sigma_{ij})\\
\end{cases}
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential = "3SPN2ExcludedVolume"
ignore.particles_within.bond       = 1
ignore.particles_within.angle      = 1
ignore.particles_within.dihedral   = 1
ignore.particles_within.nucleotide = 1
spatial_partition = {type = "CellList", margin = 0.4}
parameters = [
    {index =   0, kind = "S"},
    {index =   1, kind = "A"},
    {index =   2, kind = "P"},
    {index =   3, kind = "S"},
]
```

## Input Reference

- `index`: Integer
  - The index of the particle.
- `offset`: Integer (ptional. By default, 0.)
  - Offset value of the index.
- `kind`: String
  - One of the `"S"`(sugar), `"A", "T", "C", "G"`(base), or `"P"`(phosphate).

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
