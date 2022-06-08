+++
title = "StoichiometricUniformCubicPan"
weight = 100000
+++

# StoichiometricUniformCubicPanPotential

鍋型のポテンシャルを表現するためのシンプルな3次関数型ポテンシャルです。`StoichiometricInteraction`で使うポテンシャルです。

{{<katex display>}}
U(r) = \begin{cases}
-1 & 0 < r \leq d \\
- \left(1 - 3 \left(\frac{r - d}{\Delta r} \right)^2 + 2 \left( \frac{r - d}{\Delta r} \right)^3 \right) & d < r \leq + \Delta r\\
0 & d + \Delta r < r
\end{cases}
{{</katex>}}

全ての粒子が同じパラメータを持っている場合のためのポテンシャルです。

## 例

```toml
[[forcefields.global]]
interaction              = "Stoichiometric"
potential                = "StoichiometricUniformCubicPan"
ignore.molecule          = "Nothing"
spatial_partition.type   = {type = "CellList", margin = 0.2}
spatial_partition.margin = 0.2
particle_kinds = [
    {name = "A", coef = 1},
    {name = "B", coef = 1}
]

epsilon = 5.0  # parameter for StoichiometricInteraction
v0      = 10.0 # parameter for StoichiometricUniformCubicPanPotential
range   = 10.0 # parameter for StoichiometricUniformCubicPanPotential
paramters = [
    {index = 0, kind = "A"},
    # ...
]
```

## 入力

- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型 (省略可能)
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `v0`: 浮動小数点型
  - 鍋型の底の部分の幅を指定します。
- `range`: 浮動小数点型
  - 鍋型のふちの部分の幅を指定します。
