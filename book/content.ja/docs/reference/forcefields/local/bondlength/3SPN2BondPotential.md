+++
title = "3SPN2Bond"
weight = 1000000
+++

# 3SPN2BondPotential

3SPN2ポテンシャルの一部分です。

{{<katex display>}}
U(v) = k(v - v_0)^2 + 100k (v - v_0)^4
{{</katex>}}

## 例

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "3SPN2Bond"
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100, v0 = 1.0, k = 10.0},
    # ...
]
```

## 入力

- `k`: 浮動小数点数型
  - パラメータの強さを指定します。
- `v0`: 浮動小数点数型
  - 最安定距離を指定します。
- `indices`: 整数の配列型（長さ: 2）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。
