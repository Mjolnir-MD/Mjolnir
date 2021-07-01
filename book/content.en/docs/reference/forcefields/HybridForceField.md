+++
title = "HybridForceField"
weight = 6000
+++

# HybridForceField

HybridForceField connects two different forcefields.

{{< katex display >}}
V = \lambda V_1 + (1 - \lambda) V_2
{{< /katex >}}

The topology can be different.

## Example input

First, you need to define and name several forcefields to be connected.

Then, set `forcefields.type = "Hybrid"` in `[simulator]` table and provide `lambda` value.
The first and the second `[[forcefields]]` table will be the {{<katex>}} V_1 {{</katex>}} and {{<katex>}} V_2 {{</katex>}} in the above equation, respectively.

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "Unlimited"
precision     = "double"
delta_t       = 0.1
total_step    = 1000000
save_step     =   1_000
seed          = 2859805901
integrator.type = "BAOABLangevin"
integrator.gammas = [
    # ...
]
forcefields.type = "Hybrid"
forcefields.lambda = 0.5
# ...

[[forcefields]] # V1
[[forcefields.local]]
# ...

[[forcefields]] # V2
[[forcefields.local]]
# ...
```
