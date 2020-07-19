+++
title = "Cosine"
weight = 100000
+++

# CosinePotential

HarmonicよりもDihedralAngleで用いやすい、周期的なポテンシャルです。

{{<katex display>}}
U(v) = k\left(1 + \cos(n(v - v_0))\right)
{{</katex>}}

## 例

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Cosine"
topology    = "none"
parameters = [
    {indices = [0, 1, 2, 3], v0 = 1.57, k = 10.0, n = 1},
]
```

## 入力

- `k`: 浮動小数点数型
  - このポテンシャルの強さを指定します。
- `n`: 整数
  - 谷の個数を決めます。
- `v0`: 浮動小数点数
  - 谷の位置を決めます。
- `indices`: 整数の配列型（長さ: 4）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。

