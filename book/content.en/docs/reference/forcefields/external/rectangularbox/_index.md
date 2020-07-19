+++
title = "RectangularBox"
weight = 20000
bookCollapseSection = true
+++

# RectangularBox

`RectangularBox` interaction keeps particles inside a box.

{{<hint warning>}}
If any particle locates outside of the box at the initial configuration, it fails to start.
{{</hint>}}

{{<hint warning>}}
This interaction does not work under the periodic boundary condition.
{{</hint>}}


## Example

```toml
[[forcefields.external]]
interaction = "RectangularBox"
potential   = "ExcludedVolumeWall"

box.lower   = [  0.0,   0.0,   0.0]
box.upper   = [100.0, 100.0, 100.0]
box.margin  = 0.4

# potential related
epsilon     = 0.1
parameters  = [
    {index = 0, radius = 1.0}, # required parameters depend on potential.
]
```

## Input reference

- `interaction`: String
  - Name of the interaction. Here, it is `"RectangularBox"`.
- `potential`: String
  - [`"LennardJonesWall"`]({{<relref "LennardJonesWallPotential.md">}})
  - [`"ExcludedVolumeWall"`]({{<relref "ExcludedVolumeWallPotential.md">}})
- `box`: Table
  - `box.lower`: Array of Floats
    - lower boundary of the box.
  - `box.upper`: Array of Floats
    - upper boundary of the box.
  - `box.margin`: Float
    - margin of the neighboring list, relative to the cutoff length.
