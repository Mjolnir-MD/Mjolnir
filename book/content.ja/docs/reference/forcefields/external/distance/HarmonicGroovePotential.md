+++
title  = "HarmonicGroove"
weight = 200000
+++

# HarmonicGroovePotential

シンプルな調和振動子ポテンシャルです。

{{<katex display>}}
U(r) = k (v - v_0)^2
{{</katex>}}

## Example

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "HarmonicGroove"
shape.name     = "AxisAlignedPlane"
shape.axis     = "X"
shape.position = 0.0
shape.margin   = 1.0

parameters = [
    {index = 0, k = 1.0, v0 = 1.0},
    # ...
]
```

## Input Reference

- `k`: 浮動小数点数型
    - ポテンシャルの強さを指定します。
- `v0`: 浮動小数点数型
    - 最安定点を指定します。
