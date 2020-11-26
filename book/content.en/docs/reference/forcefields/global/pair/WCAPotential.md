+++
title = "WCA"
weight = 700000
+++

# WCAPotential

The well-known WCA (Weeks-Chandler-Andersen) potential. It is 

{{<katex display>}}
U(r) =
\begin{cases}
4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right] + \epsilon & (r < \sigma_{ij}\sqrt[6]{2})\\
0 & (r \geq \sigma_{ij}\sqrt[6]{2})\\
\end{cases}
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "WCA"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

parameters = [
    {index = 0, offset = 100, sigma = 2.0, epsilon = 10.0},
]
```

## Input Reference

To calculate {{<katex>}} \sigma {{</katex>}} and {{<katex>}} \epsilon {{</katex>}} for each pair, Lorentz-Berthelot combining rules are used.

{{<katex display>}}
\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2} \\
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j}
{{</katex>}}

- `index`: Integer
  - The index of the particle.
- `offset`: Integer (optional. By default, 0.)
  - Offset value of the index.
- `sigma`: Floating
  - It determines the effective particle size.
- `epsilon`: Floating
  - It determines the strength of the potential.

Since this potential becomes exactly 0 at {{<katex>}} r = \sigma\sqrt[6]{2} {{</katex>}}, always this cutoff distance is used. You don't need to set `cutoff`.

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
