+++
title = "Distance"
weight = 30000
bookCollapseSection = true
+++

# ExternalDistance

It depends on the distance between particles and a structure such as a planner surface.

## Example

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "LennardJonesWall"
shape.name     = "AxisAlignedPlane"
shape.axis     = "X"
shape.position = 0.0
shape.margin   = 0.5
parameters     = [
# ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"Distance"`.
- `potential`: String
  - The following potentials are available.
  - [`"LennardJonesWall"`]({{<relref "LennardJonesWallPotential.md">}})
  - [`"ExcludedVolumeWall"`]({{<relref "ExcludedVolumeWallPotential.md">}})
  - [`"HarmonicGroove"`]({{<relref "HarmonicGroovePotential.md">}})
  - [`"ImplicitMembrane"`]({{<relref "ImplicitMembranePotential.md">}})
- `shape`: Table
  - The shape and the position of a structure that interacts with particles.

### `shape`

- `name`: String
  - The following shape is available.
  - `"AxisAlignedPlane"` requires the following fields.
    - `axis`: String
      - The axis that is parpendicular to the plane. like: `"X", "Y", "Z"`.
    - `position`: Floating
      - The position of the plane along the axis specified.
- `margin`: Floating
  - The margin of the neighborling list relative to the cutoff length.
