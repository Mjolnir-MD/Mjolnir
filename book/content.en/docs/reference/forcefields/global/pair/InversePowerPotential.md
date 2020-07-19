+++
title = "InversePower"
weight = 500000
+++

# InversePowerPotential

It is a generalized excluded volume potential.

{{<katex display>}}
U(r) = \epsilon\left(\frac{\sigma}{r}\right)^{n}
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "InversePower"

ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
ignore.groups.intra = ["chain-A"]
ignore.groups.inter = [["chain-B", "chain-C"]]

spatial_partition.type = "CellList"
spatial_partition.margin = 0.2

cutoff     = 2.5
epsilon    = 0.2
n          = 5
parameters = [
    {index = 0, radius = 2.0}
]
```

## Input Reference

- `cutoff`: Floating (Optional. By default, {{<katex>}}2^{\frac{12}{n}}{{</katex>}}.)
  - Cutoff distance relative to the {{<katex>}} \sigma_{ij} {{</katex>}}.
- `epsilon`: Floating
  - It determines the strength of the potential.
- `n`: Integer
  - It determines the slope of the potential.
- `index`: Integer
  - The index of the particle. The index is 0-based.
- `offset`: Integer (Optional. By default, 0)
  - Offset for the index.
- `radius`: Floating
  - It determines the effective size of the particle.
  - {{<katex>}} \sigma {{</katex>}} in the above equation is a sum of the radii.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
