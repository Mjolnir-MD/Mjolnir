+++
title = "ExcludedVolume"
weight = 400000
+++

# ExcludedVolume

A simple excluded volume potential.

{{<katex display>}}
U(r) = \epsilon\left(\frac{\sigma}{r}\right)^{12}
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"

ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1

spatial_partition.type = "CellList"
spatial_partition.margin = 0.2

cutoff     = 2.5
epsilon    = 0.2
parameters = [
    {index = 0, radius = 2.0}
]
```

## Input Reference

- `cutoff`: Floating (Optional. By default, 2.)
  - Cutoff distance relative to the {{<katex>}} \sigma_{ij} {{</katex>}}.
- `epsilon`: Floating
  - It determines the strength of the potential.
- `index`: Integer
  - The index of the particle. The index is 0-based.
- `offset`: Integer (Optional. By default, 0)
  - Offset for the index.
- `radius`: Floating
  - It determines the effective size of the particle.
  - {{<katex>}} \sigma {{</katex>}} in the above equation is a sum of the radii.
