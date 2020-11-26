+++
title = "iSoLF"
weight = 100000
+++

# iSoLFAttractivePotential

iSoLFAttractive attractive part of the iSoLF potential for coarse-grained lipids developed by the following paper.

- Diego Ugarte La Torre and Shoji Takada (2020) J. Chem. Phys 153, 205101
  - https://doi.org/10.1063/5.0026342

{{<katex display>}}
U(r) =
\begin{cases}
-\epsilon_{ij} & r_{ij} < \sqrt[6]{2}\sigma_{ij}\\
-\epsilon_{ij} \cos^2\left[\frac{\pi}{2\omega_{ij}}(r_{ij} - \sqrt[6]{2}\sigma_{ij}) \right] & (\sqrt[6]{2}\sigma_{ij} < r_{ij} < \sqrt[6]{2}\sigma_{ij} + \omega_{ij})\\
0 & (\sqrt[6]{2}\sigma_{ij} + \omega_{ij} < r_{ij})
\end{cases}
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "iSoLFAttractive"
ignore.particles_within = {bond = 1, angle = 1}
ignore.group.inter = [
    ["T1", "T3"]
]
spatial_partition = {type = "CellList", margin = 0.5}
env.popc_epsilon = 0.416
env.popc_omega   = 9.867
env.popc_sigma_T = 7.111
parameters = [
{index =     2, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     3, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     4, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     7, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     8, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     9, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
# ...
]
```

## Input Reference

To calculate {{<katex>}} \sigma {{</katex>}}, {{<katex>}} \epsilon {{</katex>}} and {{<katex>}} \omega {{</katex>}} for each pair, Lorentz-Berthelot combining rules are used.

{{<katex display>}}
\sigma_{ij}   = \frac{\sigma_i + \sigma_j}{2} \\
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j} \\
\omega_{ij}   = \frac{\omega_i\omega_j}{2}
{{</katex>}}

- `index`: Integer
  - The index of the particle.
- `offset`: Integer (optional. By default, 0.)
  - Offset value of the index.
- `sigma`: Floating
  - It determines the effective particle size.
- `epsilon`: Floating
  - It determines the strength of the potential.
- `omega`: Floating
  - It determines the width of the attractive well.

Since this potential becomes exactly 0 at {{<katex>}} r = \sqrt[6]{2}\sigma_{ij} + \omega_{ij} {{</katex>}}, always this cutoff distance is used. You don't need to set `cutoff`.

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
