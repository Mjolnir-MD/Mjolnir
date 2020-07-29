+++
title = "LennardJones"
weight = 100000
+++

# LennardJonesWallPotential

It is a 6-12 Lennard Jones potential between particle and a spatial structure.

{{<katex display>}}
U(r) = 4\epsilon\left(\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right)
{{</katex>}}

## Example

```toml
[[forcefields.external]]
interaction = "RectangularBox"
box.lower   = [  0.0,   0.0,   0.0]
box.upper   = [100.0, 100.0, 100.0]
box.margin  = 0.4

potential   = "LennardJonesWall"
parameters  = [
    {index = 0, sigma = 1.0, epsilon = 0.1},
    # ...
]
```

## Input Reference

- `sigma`: Floating
  - It determines the effective size of the particle.
- `epsilon`: Floating
  - It determines the strength of the potential.
