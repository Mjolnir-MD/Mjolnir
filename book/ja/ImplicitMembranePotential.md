# ImplicitMembranePotential

空間中の面と粒子の間に働くImplicit Membraneポテンシャルです。

## 例

`shape`パラメータについては[ExternalDistanceInteraction](ExternalDistanceInteraction.md)を参照して下さい。

```toml
[[forcefields.external]]
interaction    = "Distance"
shape.name     = "AxisAlignedPlane"
shape.axis     = "+X"
shape.position = 0.0
shape.margin   = 0.5

potential   = "ImplicitMembrane"
thickness   = 4.0
bend        = 2.0
interaction_magnitude = 1.0
parameters  = [
    {index = 0, hydrophobicity = 0.1},
    # ...
]
```

## 入力

- `thickness`: 浮動小数点数型
  - 膜の分厚さを決めるパラメータです。粒子によらず一定です。
- `bend`: 浮動小数点数型
  - 膜の端点でのポテンシャルの傾きを調整するパラメータです。粒子によらず一定です。
- `interaction_magnitude`: 浮動小数点数型
  - ポテンシャルの強さを調整するパラメータです。粒子によらず一定です。
- `hydrophobicity`: 浮動小数点数型
  - 各粒子の疎水性のパラメータです。親水性の場合は負の値になります。

