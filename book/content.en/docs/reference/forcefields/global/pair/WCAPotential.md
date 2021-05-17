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

There are two different way to define the parameters.
You can either use the normal combining rule or define all the pair-paremeters.

By defining parameters for each particle, it uses Lorentz-Berthelot combining rule.

{{<katex display>}}
\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2} \\
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j}
{{</katex>}}

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


- `index`: Integer
  - The index of the particle.
- `offset`: Integer (optional. By default, 0.)
  - Offset value of the index.
- `sigma`: Floating
  - It determines the effective particle size. If no table is given, this is required.
- `epsilon`: Floating
  - It determines the strength of the potential. If no table is given, this is required.
- `name`: String
  - Name of the particle. If `table` is given, this is required.

Since this potential becomes exactly 0 at {{<katex>}} r = \sigma\sqrt[6]{2} {{</katex>}}, always this cutoff distance is used. You don't need to set `cutoff`.

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
