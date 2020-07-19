+++
title = "PositionRestraint"
weight = 10000
+++

# PositionRestraint

It restrains a particle to a fixed point in a space via harmonic potential.

## Example

```toml
[[forcefields.external]]
interaction = "PositionRestraint"
potential   = "Harmonic"
parameters  = [
    {index = 0, position = [0.0, 0.0, 0.0], k = 0.1, v0 = 10.0},
    # ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaciton. Here, it is `"PositionRestraint"`.
- `potential`: String
  - Name of the potential. Currently, only the harmonic potential is available.
- `parameters`: Table
  - `index`: Integer
    - The index of the particle. The index is 0-based.
  - `offset`: Integer (Optional. By default, 0.)
    - The offset value for the index.
  - `position`: Array of Floatings
    - The cartesian coordinate (x, y, and z in this order) of the position.
  - `k`: Floating
    - The strength of the restraint.
  - `v0`: Floating
    - The native distance from the point.
    - To restrain a particle exactly at the position, set 0.
