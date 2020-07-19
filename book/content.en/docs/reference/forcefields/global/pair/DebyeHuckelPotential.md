+++
title = "DebyeHückel"
weight = 300000
+++

# DebyeHückelPotential

It is a electrostatic potential in an ionic solution based on Debye-Hückel equation.

{{<katex display>}}
U(r_{ij}) = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_k r_{ij}} \exp(-r_{ij}/\lambda_D)
{{</katex>}}

{{<katex display>}}
\lambda_D = \sqrt{\frac{\epsilon_0\epsilon_k}{2\beta N_A e_c^2 I}}
{{</katex>}}

The permittivity of the water is derived from the same equation used in the following papers.

- Sambriski, E. J. et al., (2009) Biophys. J.
- Hinckley, D. M. et al., (2013) JCP.
- Freeman, G. S., et al., (2014) JCP.

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "DebyeHuckel"
ignore.particles_within.bond = 3
spatial_partition.type = "CellList"
spatial_partition.margin = 0.2
cutoff     = 5.5
parameters = [
    {index = 0, charge = 1.0},
]
```

## Input Reference

This potential requires [system attribute]({{<relref "/docs/reference/system#attribute">}}), `temperature` and `ionic_strength`.

- `cutoff`: Floating (Optional. By default, 5.5.)
  - Cutoff length relative to the debye length.
- `index`: Integer
  - The index of the particle.
- `offset`: Integer (ptional. By default, 0.)
  - Offset value of the index.
- `charge`: Floating
  - The value of charge of the particle.
