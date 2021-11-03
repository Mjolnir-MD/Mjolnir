+++
title = "CoMPullingForce"
weight = 60000
+++

# CoMPullingForce

It applies a force to the center of mass of the specified particles.

## Example

```toml
[[forcefields.external]]
interaction = "CoMPullingForce"
parameters  = [
    {indices = [0, 1, 2], force = [0.0, 0.0, 0.0144]},
    {indices = "[0, 10)", force = 0.0144, direction = [1.0, 1.0, 1.0]},
    {indices = ["[0, 99]", "[200, 300)"], force = 0.0144, direction = [1.0, 1.0, 1.0]},
    # ...
]
```

## Input Reference

- `interaction`: String
  - Name of the interaciton. Here, it is `"CoMPullingForce"`.
- `parameters`: Table
  - `indices`: Array of Integers, String, or Array of Strings
    - The indices of the particles.
    - If it is an Array of Integers, then the integers will be the indices of the particles.
    - You can indicate an interval of the indices by a standard form. A string "[0, 10)" means integers from 0 to 9 (half-open interval).
    - If it is an Array of Strings, then the intervals will be concatenated.
    - Note that the indices will be made unique.
  - `force`: Array of Floatings (length = 3) or Floating
    - If it is an Array of Floatings, it represents the force to apply.
    - If it is Floating, it represents the strength of the force. In that case `direction` must be provided.
    - The unit depends on the unit system you choose in [`[units]`]({{<relref "../../units/_index.md" >}}) (e.g. `kcal/mol/Å`).
       - `1 kcal/mol/Å ~ 69.5 pN`.
  - `direction`: Array of Floatings
    - The cartesian coordinate (x, y, and z in this order) of the force direction. It will be normalized.
