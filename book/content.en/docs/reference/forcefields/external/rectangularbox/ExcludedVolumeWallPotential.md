+++
title = "ExcludedVolume"
weight = 200000
+++

# ExcludedVolumeWallPotential

It is a excluded volume potential between particle and a spatial structure.

{{<katex display>}}
U(r) = \epsilon\left(\frac{r_0}{r}\right)^{12}
{{</katex>}}

## Example

```toml
[[forcefields.external]]
interaction = "RectangularBox"
box.lower   = [  0.0,   0.0,   0.0]
box.upper   = [100.0, 100.0, 100.0]
box.margin  = 0.4

potential   = "LennardJonesWall"
epsilon     = 0.5
parameters  = [
    {index = 0, radius = 1.0},
    # ...
]
```

## Input Reference

- `epsilon`: Floating
  - It determines the strength of the potential.
- `radius`: Floating
  - It determines the effective radius of the particle.
