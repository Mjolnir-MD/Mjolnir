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
interaction = "RectangularBox"
box.lower   = [  0.0,   0.0,   0.0]
box.upper   = [100.0, 100.0, 100.0]
box.margin  = 0.4

potential   = "LennardJonesWall"
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
