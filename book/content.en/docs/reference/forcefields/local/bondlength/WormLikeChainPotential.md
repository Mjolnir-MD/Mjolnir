+++
title = "WormLikeChain"
weight = 600000
+++

# WormLikeChain

A potential function based on the Worm-Like chain model.

{{<katex display>}}
U(r) = \frac{k_B T}{p}  \left(\frac{l_c}{4} \left[ \frac{1}{1 - \frac{r}{l_c}} - 1 \right] - \frac{r}{4} + \frac{r^2}{2l_c}\right)
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "WormLikeChain"
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100, p = 5.0, lc = 100.0},
    # ...
]
```

## Input Reference

- `p`: Floating
  - The persistent length of the polymer.
- `lc`: Floating
  - The maximum length of the polymer.
- `indices`: Array of Integers (length = 2)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
