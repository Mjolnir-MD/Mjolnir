# ConstraintForceField



## Example

```toml
[[forcefields.constraint]]
topology      = "bond"
max_iteration = 500
tolerance     = 1e-6
parameters  = [
    {indices = [0, 1], offset = 100, v0 = 3.8},
]
```

## Input reference

- `topology`: String
  - Name of this connection in [`"Topology"`](Topology.md). It affects to the global interactions.
- `max_iteration`: Integer
  - The maximum number of iteration while correcting the positions of particles.
  - If it does not converge after the iteration, it outputs warning.
- `tolerance`: Float
  - The acceptable error range in the distance.
- `parameters`: Array of Tables
  - `indices`: Array of Integers
    - The indices of the particles that are constrained. The indices are 0-based.
  - `offset`: Integer(optional)
    - Offset of the index. It can be helpful if there are some groups of particles.
  - `v0`: Float
    - The native distance.
