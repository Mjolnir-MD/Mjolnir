+++
title = "ImplicitMembrane"
weight = 300000
+++

# ImplicitMembranePotential

平面の周りにある粒子をその疎水度に応じて安定化するimplicitな膜モデルです。

{{<katex display>}}
U(\mathbf{r}) = \sum_i^N k h_i \tanh\left(\mathrm{bend} * \left(|z_i - z_0| - \frac{\mathrm{thickness}}{2}\right)\right)
{{</katex>}}

## 例

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "ImplicitMembrane"
shape.name     = "AxisAlignedPlane"
shape.axis     = "X"
shape.position = 0.0
shape.margin   = 0.5

bend           = 1.0
thickness      = 4.0
interaction_magnitude = 10.0
parameters  = [
    {index = 0, hydrophobicity = 1.0},
    # ...
]
```

## 入力

- `bend`: 浮動小数点数
  - 中央の安定化を受ける部分と外側の部分との間のスロープの傾きを決めます。
- `thickness`: 浮動小数点数
  - 中央の安定化を受ける部分の幅を決めます。
- `interaction_magnitude`: 浮動小数点数
  - ポテンシャルの強さを決めます。上式の{{<katex>}} k {{</katex>}}です。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 粒子の番号です。最初の粒子は0番目です。
  - `hydrophobicity`: 浮動小数点数型
    - 粒子の疎水度です。上式の{{<katex>}} h {{</katex>}}です。

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
