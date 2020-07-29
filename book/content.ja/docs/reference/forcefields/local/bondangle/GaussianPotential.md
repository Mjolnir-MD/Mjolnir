+++
title = "Gaussian"
weight = 200000
+++

# GaussianPotential

ガウシアンポテンシャルは以下のような形のポテンシャルです。

{{<katex display>}}
U(v) = k\exp\left(\frac{-(v - v_0)^2}{2\sigma^2}\right)
{{</katex>}}

## 例

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Gaussian"
topology    = "none"
parameters  = [
    {indices = [0, 1, 2], v0 = 1.0, k = -100.0, sigma = 5.0},
    # ...
]
```

## 入力

- `k`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `sigma`: 浮動小数点数型
  - ポテンシャルが効果を及ぼす幅を指定します。
- `v0`: 浮動小数点数型
  - 最安定点を指定します。
- `indices`: 整数の配列型（長さ:3）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。
