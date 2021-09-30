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

## make λ dynamic

We can integrate {{< katex >}}\lambda{{< /katex >}} over time by considering it as a dynamic variable.
To do that, set `forcefields.lambda = "dynamic"` and define `dynamic_variables.lambda` in `[[systems]]`.
The coordinate, velocity, force are saved in the `.ene` file.

Also, you can apply a harmonic potential {{< katex >}}k(\lambda - \lambda_0)^2{{< /katex >}} to {{< katex >}}\lambda{{< /katex >}}.

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
forcefields.lambda = "dynamic" # make lambda dynamic
forcefields.v0_lambda = 1.0    # k(v - v0)^2
forcefields.k_lambda  = 100.0  # k(v - v0)^2

# ...

[[systems]]
attributes.temperature = 300.0
dynamic_variables.lambda = {x = 1.0, m = 1000.0, gamma = 1e-5, boundary = "Repulsive", lower = 0.0, upper = 1.0}

[[forcefields]] # V1
[[forcefields.local]]
# ...

[[forcefields]] # V2
[[forcefields.local]]
# ...
```
