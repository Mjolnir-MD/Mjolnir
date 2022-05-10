+++
title = "LennardJonesAttractive"
weight = 800000
+++

# LennardJonesAttractivePotential

The attractive part of the Lennard-Jones potential.

{{<katex display>}}
U(r) = \begin{cases}
-\epsilon & (r < \sigma_{ij}\sqrt[6]{2})\\
4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right] & (\sigma_{ij}\sqrt[6]{2} < r)
\end{cases}
{{</katex>}}

## Example

There are two different way to define the parameters.
You can either use the normal combining rule or define all the pair-paremeters.

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJonesAttractive"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

cutoff = 2.5
parameters = [
    {index = 0, offset = 100, sigma = 2.0, epsilon = 10.0},
    # ...
]
```

To provide pair-parameters manually, define `table` and give `name`s to the particles.

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJonesAttractive"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}
cutoff = 2.5

table.A.A = {sigma = 1.0, epsilon = 2.0}
table.A.B = {sigma = 3.0, epsilon = 1.0} # B.A will be the same
table.B.B = {sigma = 1.0, epsilon = 1.5}
parameters = [
    {index = 0, offset = 100, name = "A"},
    {index = 1, offset = 100, name = "B"},
    # ...
]
```

## Input Reference

- `cutoff`: Floating (Optional. By default, 2.5.)
  - The cutoff distance relative to the maximum {{<katex>}}\sigma_{ij}{{</katex>}}.
- `index`: Integer
  - The index of the particle.
- `offset`: Integer (ptional. By default, 0.)
  - Offset value of the index.
- `sigma`: Floating
  - It determines the effective particle size. If no table is given, this is required.
- `epsilon`: Integer
  - It determines the strength of the potential. If no table is given, this is required.
- `name`: String
  - Name of the particle. If `table` is given, this is required.

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
