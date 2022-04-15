+++
title  = "UniformCubicPan"
weight = 200000
+++

# UniformCubicPanPotential

鍋型のポテンシャルを表現するためのシンプルな3次関数型ポテンシャルです。

{{<katex display>}}
U(r) = \begin{cases}
-\varepsilon & 0 < r \leq d \\
-\varepsilon \left(1 - 3 \left(\frac{r - d}{\Delta r} \right)^2 + 2 \left( \frac{r - d}{\Delta r} \right)^3\right) & d < r \leq d + \Delta r\\
0 & d + \Delta r < r
\end{cases}
{{</katex>}}

全ての粒子が同じパラメータを持っている場合のためのポテンシャルです。

## 例

```
[[forcefields.global]]
interaction = "Pair"
potential   = "CubicPanPotential"
ignore.molecule = "Nothing"
ignore.partcles_within.bond     = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

epsilon = 1.0
v0      = 9.0
range   = 10.0
parameters = [
    {index = 0, offset = 100}, # to control which particle participantes
]
```

## 入力

- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型 (省略可能)
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `epsilon`: 浮動小数点型
  - ポテンシャルの強さを指定します。
- `v0`: 浮動小数点型
  - 鍋型の底の部分の幅を指定します。
- `range`: 浮動小数点型
  - 鍋型のふちの部分の幅を指定します。
