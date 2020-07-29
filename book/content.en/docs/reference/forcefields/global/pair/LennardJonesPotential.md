+++
title = "LennardJones"
weight = 100000
+++

# LennardJonesPotential

The well-known Lennard-Jones potential.

{{<katex display>}}
U(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
{{</katex>}}

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

cutoff = 2.5
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

- `cutoff`: Floating (Optional. By default, 2.5.)
  - The cutoff distance relative to the maximum {{<katex>}}\sigma_{ij}{{</katex>}}.
- `index`: Integer
  - The index of the particle.
- `offset`: Integer (ptional. By default, 0.)
  - Offset value of the index.
- `sigma`: Floating
  - It determines the effective particle size.
- `epsilon`: Integer
  - It determines the strength of the potential.

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
