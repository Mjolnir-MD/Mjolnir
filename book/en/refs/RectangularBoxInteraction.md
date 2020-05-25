# RectangularBox

`RectangularBox` interaction keeps particles inside a box.

{% hint style='info' %}
If any particle locates outside of the box at the initial configuration, it fails to start.
{% endhint %}

{% hint style='info' %}
This interaction does not work under the periodic boundary condition.
{% endhint %}


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
    {index = 0, radius = 1.0},
]
```

## Input reference

- `interaction`: String
  - `"RectangularBox"`.
- `potential`: String
  - [`"LennardJonesWall"`](LennardJonesWallPotential.md)
  - [`"ExcludedVolumeWall"`](ExcludedVolumeWallPotential.md)
- `box`: Table
  - `box.lower`: Array of Floats
    - lower boundary of the box.
  - `box.upper`: Array of Floats
    - upper boundary of the box.
  - `box.margin`: Float
    - margin of the neighboring list

Depending on the potentials, required paramters may change.
See potential pages for detail.
