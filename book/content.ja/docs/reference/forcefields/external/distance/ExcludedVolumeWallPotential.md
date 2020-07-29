+++
title = "ExcludedVolume"
weight = 200000
+++

# ExcludedVolumeWallPotential

排除体積相互作用のポテンシャルです。

{{<katex display>}}
U(r) = \epsilon\left(\frac{r_0}{r}\right)^{12}
{{</katex>}}

## 例

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "ExcludedVolumeWall"
shape.name     = "AxisAlignedPlane"
shape.axis     = "X"
shape.position = 0.0
shape.margin   = 0.5

epsilon     = 0.5
parameters  = [
    {index = 0, radius = 1.0},
    # ...
]
```

## 入力

- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを決めます。
- `radius`: 浮動小数点数型
  - 粒子の半径を決めます。
