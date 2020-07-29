+++
title = "FlexibleLocal"
weight = 300000
+++

# FlexibleLocalAngle

It represents angle distribution of flexible regions of proteins in a Coarse-Grained protein model.
Since it represents a distribution, it is based on values on grid points and interpolates those values.

It is a component of AICG2+ Coarse-Grained protein model.

It is developed in the following paper.

- T. Terakawa and S. Takada, (2011) Biophys J.

## Example

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential = "FlexibleLocalAngle"
topology  = "none"
env.x      = [ 1.3090,  1.4835,  1.6581,  1.8326,  2.0071,  2.1817,  2.3562,  2.5307,  2.7053,  2.8798]
env.y1_ALA = [   5.00,    1.34,    0.84,    1.17,    0.82,    1.00,    1.27,    1.52,    3.20,   10.00]
env.y2_ALA = [   0.00,  151.96,   14.61,  -46.89,   39.04,   -4.86,   -1.86,    8.38,  250.03,    0.00]
parameters = [
    {indices = [0, 1, 2], k = 1.0, x = "x", y = "y1_ALA", d2y = "y2_ALA"},
    # ...
]
```

## Input Reference

- `k`: Floating
  - It determines the strength of the potential.
- `x`: Array of Floating (length = 10. Optional. By default, use default value in AICG2+).
  - Grid points on which the values are specified. The unit is the radian.
  - It should be equidistant.
- `y`: Array of Floating (length = 10)
  - Potential energy value at each grid point.
- `d2y`: Array of Floating (length = 10)
  - Second derivative of potential energy value at each grid point.
- `indices`: Array of Integers (length = 3)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.
