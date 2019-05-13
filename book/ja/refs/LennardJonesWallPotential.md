# LennardJonesWallPotential

空間中の形状と粒子の間に適用される、(6,12)型のレナード・ジョーンズポテンシャルです。

## 例

`shape`パラメータについては[ExternalDistanceInteraction](ExternalDistanceInteraction.md)を参照して下さい。

```toml
[[forcefields.external]]
interaction    = "Distance"
shape.name     = "AxisAlignedPlane"
shape.axis     = "+X"
shape.position = 0.0
shape.margin   = 0.5

potential   = "LennardJonesWall"
parameters  = [
    {index = 0, sigma = 1.0, epsilon = 0.1},
    # ...
]
```

## 入力

- `sigma`: 浮動小数点数型
  - 粒子径に相当するパラメータです。
- `epsilon`: 浮動小数点数型
  - エネルギーに相当するパラメータです。

