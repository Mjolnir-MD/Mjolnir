+++
title = "PullingForce"
weight = 50000
+++

# PullingForce

It applies a force to the specified particle.

## Example

```toml
[[forcefields.external]]
interaction = "PullingForce"
parameters  = [
    {index = 0, force = [1.0, 0.0, 0.0]},
    {index = 0, force = 1.0, direction = [1.0, 1.0, 1.0]},
    # ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaciton. Here, it is `"PullingForce"`.
- `parameters`: Table
  - `index`: Integer
    - The index of the particle. The index is 0-based.
  - `force`: Array of Floatings (length = 3) or Floating
    - If it is an Array of Floatings, it represents the force to apply.
    - If it is Floating, it represents the strength of the force. In that case `direction` must be provided.
    - The unit depends on the unit system you choose in [`[units]`]({{<relref "../../units/_index.md" >}}) (e.g. `kcal/mol/Å`).
       - `1 kcal/mol/Å ~ 69.5 pN`.
  - `direction`: Array of Floatings
    - The cartesian coordinate (x, y, and z in this order) of the force direction. It will be normalized.
