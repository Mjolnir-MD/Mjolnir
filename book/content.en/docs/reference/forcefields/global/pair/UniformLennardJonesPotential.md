+++
title = "UniformLennardJones"
weight = 200000
+++

# UniformLennardJonesPotential

The well-known Lennard-Jones potential.

{{<katex display>}}
U(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
{{</katex>}}

It is a special case where all the particles has the same parameter.

## Example

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "UniformLennardJones"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

cutoff  = 2.5
sigma   = 2.0
epsilon = 0.5
parameters = [
    {index = 0, offset = 100}, # to control which particle participates
]
```

## Input Reference

- `cutoff`: Floating (Optional. By default, 2.5.)
  - The cutoff distance relative to the maximum {{<katex>}}\sigma_{ij}{{</katex>}}.
- `index`: Integer
  - The index of the particle.
- `offset`: Integer (ptional. By default, 0.)
  - Offset of the index.
- `sigma`: Floating
  - It determines the effective particle size.
- `epsilon`: Integer
  - It determines the strength of the potential.

For other values, see [Pair]({{<relref "/docs/reference/forcefields/global/pair">}}).
