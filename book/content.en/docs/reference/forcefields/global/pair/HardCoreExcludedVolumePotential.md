+++
title = "HardCore"
weight = 600000
+++

# HardCoreExcludedVolume

An excluded volume potential with a hard core that never overlaps to each other.

{{<katex display>}}
U(r) = \epsilon\left(\frac{\sigma}{r - r_0}\right)^{12}
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "HardCoreExcludedVolume"

ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1

spatial_partition.type = "CellList"
spatial_partition.margin = 0.2

cutoff     = 2.5
epsilon    = 0.2
parameters = [
    {index = 0, core_radius = 3.0, soft_shell_thickness = 2.0},
    {index = 1, core_radius = 3.0, soft_shell_thickness = 2.0},
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
- `core_radius`: Floating
  - The size of the core.
  - {{<katex>}} r_0 {{</katex>}} is a sum of `core_radii`.
- `soft_shell_thickness`: Floating
  - The thickness of the soft shell.
  - {{<katex>}}\sigma_{ij}{{</katex>}} is a sum of `soft_shell_thickness`es.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
