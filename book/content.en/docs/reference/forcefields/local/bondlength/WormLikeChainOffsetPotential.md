+++
title = "WormLikeChainOffset"
weight = 600000
+++

# WormLikeChainOffset

A potential function based on the Worm-Like chain model with distance offset.

{{<katex display>}}
U(r) = \begin{cases}
0 & ( r < l_0 ) \\
\frac{k_B T}{p}  \left(\frac{l_c}{4} \left[ \frac{1}{1 - \frac{r - l_0}{l_c}} - 1 \right] - \frac{r - l_0}{4} + \frac{\left(r - l_0 \right)^2}{2l_c}\right) & ( r \geq l_0 ) \\
\end{cases}
{{</katex>}}

## Example

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "WormLikeChainOffset"
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100, p = 5.0, lc = 100.0, l0 = 30.0},
    # ...
]
```

## Input Reference

- `p`: Floating
  - The persistent length of the polymer.
- `lc`: Floating
  - The maximum length of the polymer.
- `l0`: Floating
  - The offset length of the distance between two particles.
- `indices`: Array of Integers (length = 2)
  - Indices of particles that interact with each other. The index is 0-based.
- `offset`: Integer(Optional. By default, 0.)
  - Offset of index.

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
