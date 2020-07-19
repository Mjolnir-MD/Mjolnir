+++
title = "Dummy"
weight = 60000
+++

# DummyInteraction

`DummyInteraction` does not calculate anything. No potential is available.

It is provided to add a connection between particles in its Topology without make the particles interact.
That means that, this interaction enables to ignore some of the global interactions without having local interaction.

For more detail about topology, see [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).

## Example

```toml
[[forcefields.local]]
interaction = "Dummy"
# No potential field required.
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100}, # No other parameters.
    # ...
]
```

## Input reference

- `interaction`: String
  - Name of the interaction. Here, it is `"Dummy"`.
- `topology`: String
  - Name of the connection in `Topology`.
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 2)
    - The indices of particles to be connected. The index is 0-based.
  - `offset`: Integer (Optional. By default, 0.)
    - Offset of index.
