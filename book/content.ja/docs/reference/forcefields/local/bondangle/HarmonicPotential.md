+++
title = "Harmonic"
weight = 100000
+++

# HarmonicPotential

シンプルで、様々な用途で用いられる調和振動子ポテンシャルです。

{{<katex display>}}
U(v) = k(v-v_0)^2
{{</katex>}}

## 例

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Harmonic"
topology    = "none"
parameters  = [
    {indices = [0, 1, 2], offset = 100, v0 = 1.0, k = 100.0},
    # ...
]
```

## 入力

- `k`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `v0`: 浮動小数点数型
  - 最安定点を指定します。
- `indices`: 整数の配列型（長さ:3）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。
