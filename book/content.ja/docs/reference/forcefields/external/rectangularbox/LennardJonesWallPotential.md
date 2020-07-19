+++
title = "LennardJones"
weight = 100000
+++

# LennardJonesWallPotential

よく知られたLennard-Jonesポテンシャルです。

{{<katex display>}}
U(r) = 4\epsilon\left(\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right)
{{</katex>}}

## 例

```toml
[[forcefields.external]]
interaction = "RectangularBox"
box.lower   = [  0.0,   0.0,   0.0]
box.upper   = [100.0, 100.0, 100.0]
box.margin  = 0.4

potential   = "LennardJonesWall"
parameters  = [
    {index = 0, sigma = 1.0, epsilon = 0.1},
    # ...
]
```

## 入力

- `epsilon`: 浮動小数点数
  - ポテンシャルの強さを決めます。
- `sigma`: 浮動小数点数
  - 粒子のサイズを決めます。

